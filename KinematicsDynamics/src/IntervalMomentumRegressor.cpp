#include "IntervalMomentumRegressor.h"

namespace RAPTOR {

IntervalMomentumRegressor::IntervalMomentumRegressor(const Model& model_input, 
                                                     const std::shared_ptr<TrajectoryData>& trajPtr_input,
                                                     const SensorNoiseInfo sensor_noise_input,
                                                     Eigen::VectorXi jtype_input) :
    jtype(jtype_input),
    sensor_noise(sensor_noise_input) {
    trajPtr_ = trajPtr_input;
    N = trajPtr_->N;
    modelPtr_ = std::make_shared<Model>(model_input);
    dataPtr_ = std::make_shared<Data>(model_input);

    if (jtype.size() > 0) {
        if (modelPtr_->nv != jtype.size()) {
            std::cerr << "modelPtr_->nv = " << modelPtr_->nv << std::endl;
            std::cerr << "jtype.size() = " << jtype.size() << std::endl;
            throw std::invalid_argument("modelPtr_->nv != jtype.size()");
        }
    }
    else {
        jtype = convertPinocchioJointType(*modelPtr_);
    }

    NB = modelPtr_->nv;

    Xtree.resize(1, modelPtr_->nv);

    numParams = 10 * modelPtr_->nv;

    phi.resize(numParams);

    for (int i = 0; i < modelPtr_->nv; i++) {
        const int pinocchio_joint_id = i + 1; // the first joint in pinocchio is the root joint

        // plux in Roy Featherstone's code (transformation matrix from parent to child)
        Xtree(i) = Utils::plux(modelPtr_->jointPlacements[pinocchio_joint_id].rotation().transpose(), 
                               modelPtr_->jointPlacements[pinocchio_joint_id].translation());

        phi.segment<10>(10 * i) = 
            modelPtr_->inertias[pinocchio_joint_id]
                .toDynamicParameters();
    }

    a_grav << modelPtr_->gravity.angular(),
              modelPtr_->gravity.linear();

    tau.resize(1, N);
    ptau_pz.resize(1, N);

    for (int i = 0; i < N; i++) {
        tau(i) = VecXInt::Zero(modelPtr_->nv);
        ptau_pz(i) = MatXInt::Zero(modelPtr_->nv, trajPtr_->varLength);
    }

    Y = MatXInt::Zero(N * modelPtr_->nv, numParams);

    // In the original RegressorInverseDynamics, z represents the decision variable, 
    // such as coefficients of the polynomial trajectory.
    // Here, z represents the q, q_d, q_dd themselves since we are using TrajectoryData class.
    // trajPtr_->varLength should just be 3 * modelPtr_->nv = 3 * trajPtr_->Nact.
    pY_pz.resize(trajPtr_->varLength);
    for (int i = 0; i < trajPtr_->varLength; i++) {
        pY_pz(i) = MatXInt::Zero(N * modelPtr_->nv, numParams);
    }

    Y_CTv.resize(N * modelPtr_->nv, numParams);
    pY_CTv_pz.resize(trajPtr_->varLength);
    for (int i = 0; i < trajPtr_->varLength; i++) {
        pY_CTv_pz(i).resize(N * modelPtr_->nv, numParams);
    }
};

void IntervalMomentumRegressor::compute(const VecXd& z,
                                        bool compute_derivatives,
                                        bool compute_hessian) {
    if (trajPtr_ == nullptr) {
        throw std::runtime_error("trajPtr_ is not defined yet!");
    }

    trajPtr_->compute(z, compute_derivatives);

    int i = 0;
    #pragma omp parallel for shared(trajPtr_, modelPtr_, jtype, Xtree, Y, pY_pz, tau, ptau_pz) private(i) schedule(dynamic, 1)
    for (i = 0; i < N; i++) {
        const VecXd& q_meas = trajPtr_->q(i);
        const VecXd& q_d_meas = trajPtr_->q_d(i);

        // consider sensor noise here
        VecXInt q(q_meas.size());
        VecXInt q_d(q_d_meas.size());
        for (int j = 0; j < q_meas.size(); j++) {
            q(j) = q_meas(j) + sensor_noise.position_error;
            q_d(j) = q_d_meas(j) + sensor_noise.velocity_error;
        }

        // below is the extended regressor algorithm from ROAM-Lab
        // refer to https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/dynamics/Hqd_and_CTqd.m

        Vec6Int vJ;
        Mat6Int XJ, dXJdq;
        MatRegressorInt bodyRegressor;
        Vec6Int Sd;
        MatXInt pSd_pz(6, trajPtr_->varLength);

        Eigen::Array<Mat6Int, 1, Eigen::Dynamic> Xup(modelPtr_->nv);
        Eigen::Array<Mat6Int, 1, Eigen::Dynamic> dXupdq(modelPtr_->nv);
        Eigen::Array<Vec6Int, 1, Eigen::Dynamic> S(modelPtr_->nv);
        Eigen::Array<Vec6Int, 1, Eigen::Dynamic> v(modelPtr_->nv);
        Eigen::Array<MatXInt, 1, Eigen::Dynamic> pv_pz(modelPtr_->nv);
        Eigen::Array<MatRegressorInt, 1, Eigen::Dynamic> pbodyRegressor_pz(trajPtr_->varLength);

        // forward pass
        for (int j = 0; j < modelPtr_->nv; j++) {
            const int pinocchio_joint_id = j + 1; // the first joint in pinocchio is the root joint
            const int parent_id = modelPtr_->parents[pinocchio_joint_id] - 1;

            if (compute_derivatives) {
                Spatial::jcalc(XJ, dXJdq, S(j), jtype(j), q(j));
            }
            else {
                Spatial::jcalc(XJ, S(j), jtype(j), q(j));
            }

            vJ = S(j) * q_d(j);
            Xup(j) = intervalDouble66MatrixMultiply(XJ, Xtree(j));

            if (compute_derivatives) {
                dXupdq(j) = intervalDouble66MatrixMultiply(dXJdq, Xtree(j));
            }

            if (parent_id > -1) {
                v(j) = intervalMatrixMultiply(Xup(j), v(parent_id)) + vJ;

                if (compute_derivatives) {
                    // compute pv_pz
                    pv_pz(j) = intervalMatrixMultiply(Xup(j), pv_pz(parent_id));

                    if (j < trajPtr_->Nact) {
                        // deal with S(j) * q_d(j);
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k).lower() != 0 && S(j)(k).upper() != 0) {
                                pv_pz(j)(k, j + modelPtr_->nv) += S(j)(k);
                            }
                        }

                        // deal with dXupdq(j) * v(parent_id)
                        auto dXupdq_v_parent_id = intervalMatrixMultiply(dXupdq(j), v(parent_id));
                        for (int k = 0; k < 6; k++) {
                            pv_pz(j)(k, j) += dXupdq_v_parent_id(k);
                        }
                    }
                }
            }
            else {
                v(j) = vJ;

                if (compute_derivatives) {// compute pv_pz
                    pv_pz(j) = MatXInt::Zero(6, trajPtr_->varLength);

                    if (j < trajPtr_->Nact) {
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k).lower() != 0 && S(j)(k).upper() != 0) {
                                pv_pz(j)(k, j + modelPtr_->nv) += S(j)(k);
                            }
                        }
                    }
                }
            }

            // compute local bodyRegressor
            const Interval& v1 = v(j)(0);
            const Interval& v2 = v(j)(1);
            const Interval& v3 = v(j)(2);
            const Interval& v4 = v(j)(3);
            const Interval& v5 = v(j)(4);
            const Interval& v6 = v(j)(5);

            bodyRegressor << 
                 0,     0,   v6,   -v5,    v1,    v2,     0,   v3,    0,    0,
                 0,   -v6,    0,    v4,     0,    v1,    v2,    0,   v3,    0,
                 0,    v5,  -v4,     0,     0,     0,     0,   v1,   v2,   v3,
                v4,     0,  -v3,    v2,     0,     0,     0,    0,    0,    0,
                v5,    v3,    0,   -v1,     0,     0,     0,    0,    0,    0,
                v6,   -v2,   v1,     0,     0,     0,     0,    0,    0,    0;
                 
            if (compute_derivatives) {
                for (int k = 0; k < trajPtr_->varLength; k++) {
                    const Interval& pv1 = pv_pz(j)(0, k);
                    const Interval& pv2 = pv_pz(j)(1, k);
                    const Interval& pv3 = pv_pz(j)(2, k);
                    const Interval& pv4 = pv_pz(j)(3, k);
                    const Interval& pv5 = pv_pz(j)(4, k);
                    const Interval& pv6 = pv_pz(j)(5, k);

                    pbodyRegressor_pz(k) <<
                        0,       0,   pv6,   -pv5,   pv1,   pv2,     0,  pv3,    0,    0,
                        0,    -pv6,     0,    pv4,     0,   pv1,   pv2,    0,  pv3,    0,
                        0,     pv5,  -pv4,      0,     0,     0,     0,  pv1,  pv2,  pv3,
                        pv4,     0,  -pv3,    pv2,     0,     0,     0,    0,    0,    0,
                        pv5,   pv3,     0,   -pv1,     0,     0,     0,    0,    0,    0,
                        pv6,  -pv2,   pv1,      0,     0,     0,     0,    0,    0,    0;
                }
            }

            // reverse substitution
            for (int h = j; h >= 0; h = modelPtr_->parents[h + 1] - 1) {
                // regressor for the system momentum
                Y.block(i * modelPtr_->nv + h, j * 10, 1, 10) = 
                    intervalMatrixMultiply(S(j).transpose(), bodyRegressor);

                if (compute_derivatives) {
                    for (int k = 0; k < trajPtr_->varLength; k++) {
                        pY_pz(k).block(i * modelPtr_->nv + h, j * 10, 1, 10) = 
                            intervalMatrixMultiply(S(j).transpose(), pbodyRegressor_pz(k));
                    }
                }

                Sd = intervalMatrixMultiply(Spatial::crm(v(h)), S(h));

                if (compute_derivatives) {
                    pSd_pz.row(0) = S(h)(2) * pv_pz(h).row(1) - S(h)(1) * pv_pz(h).row(2);
                    pSd_pz.row(1) = S(h)(0) * pv_pz(h).row(2) - S(h)(2) * pv_pz(h).row(0);
                    pSd_pz.row(2) = S(h)(1) * pv_pz(h).row(0) - S(h)(0) * pv_pz(h).row(1);
                    pSd_pz.row(3) = S(h)(2) * pv_pz(h).row(4) - S(h)(1) * pv_pz(h).row(5) - S(h)(4) * pv_pz(h).row(2) + S(h)(5) * pv_pz(h).row(1);
                    pSd_pz.row(4) = S(h)(0) * pv_pz(h).row(5) - S(h)(2) * pv_pz(h).row(3) + S(h)(3) * pv_pz(h).row(2) - S(h)(5) * pv_pz(h).row(0);
                    pSd_pz.row(5) = S(h)(1) * pv_pz(h).row(3) - S(h)(0) * pv_pz(h).row(4) - S(h)(3) * pv_pz(h).row(1) + S(h)(4) * pv_pz(h).row(0);
                }

                // regressor for C^T(q, v) * v
                Y_CTv.block(i * modelPtr_->nv + h, j * 10, 1, 10) = 
                    intervalMatrixMultiply(Sd.transpose(), bodyRegressor);

                if (compute_derivatives) {
                    for (int k = 0; k < trajPtr_->varLength; k++) {
                        pY_CTv_pz(k).block(i * modelPtr_->nv + h, j * 10, 1, 10) = 
                            intervalMatrixMultiply(pSd_pz.col(k).transpose(), bodyRegressor) +
                            intervalMatrixMultiply(Sd.transpose(), pbodyRegressor_pz(k));
                    }
                }

                if (modelPtr_->parents[h + 1] > 0) {
                    if (compute_derivatives) {
                        for (int k = 0; k < trajPtr_->varLength; k++) {
                            pbodyRegressor_pz(k) = 
                                intervalMatrixMultiply(Xup(h).transpose(), pbodyRegressor_pz(k));
                        }
                        pbodyRegressor_pz(h) += intervalMatrixMultiply(dXupdq(h).transpose(), bodyRegressor);
                    }

                    bodyRegressor = intervalMatrixMultiply(Xup(h).transpose(), bodyRegressor);
                }
            }
        }

        // compute torque
        tau(i) = intervalDoubleMatrixMultiply(Y.middleRows(i * modelPtr_->nv, modelPtr_->nv), phi);

        if (compute_derivatives) {
            for (int k = 0; k < trajPtr_->varLength; k++) {
                ptau_pz(i).col(k) = intervalDoubleMatrixMultiply(pY_pz(k).middleRows(i * modelPtr_->nv, modelPtr_->nv), phi);
            }
        }
    }
}

}; // namespace RAPTOR