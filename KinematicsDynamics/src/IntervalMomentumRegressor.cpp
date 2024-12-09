#include "IntervalMomentumRegressor.h"

namespace RAPTOR {

IntervalMomentumRegressor::IntervalMomentumRegressor(const Model& model_input, 
                                                     const std::shared_ptr<TrajectoryData>& trajPtr_input,
                                                     Eigen::VectorXi jtype_input) :
    jtype(jtype_input) {
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
        ptau_pz(i) = MatXInt::Zero(modelPtr_->nv, 2 * modelPtr_->nv);
    }

    Y = MatXInt::Zero(N * modelPtr_->nv, numParams);

    // In the original RegressorInverseDynamics, z represents the decision variable, 
    // such as coefficients of the polynomial trajectory.
    // Here, z represents the q, q_d themselves since we are using TrajectoryData class.
    pY_pz.resize(2 * modelPtr_->nv);
    for (int i = 0; i < 2 * modelPtr_->nv; i++) {
        pY_pz(i) = MatXInt::Zero(N * modelPtr_->nv, numParams);
    }

    Y_CTv.resize(N * modelPtr_->nv, numParams);
    pY_CTv_pz.resize(2 * modelPtr_->nv);
    for (int i = 0; i < 2 * modelPtr_->nv; i++) {
        pY_CTv_pz(i).resize(N * modelPtr_->nv, numParams);
    }

    mrPtr_ = std::make_shared<MomentumRegressor>(model_input, trajPtr_, jtype_input);
};

void IntervalMomentumRegressor::compute(const VecXd& z,
                                        bool compute_derivatives,
                                        bool compute_hessian) {
    if (trajPtr_ == nullptr) {
        throw std::runtime_error("trajPtr_ is not defined yet!");
    }

    trajPtr_->compute(z, compute_derivatives);

    const auto& sensor_noise = trajPtr_->sensor_noise;

    int i = 0;
    #pragma omp parallel for shared(trajPtr_, modelPtr_, jtype, Xtree, Y, pY_pz, tau, ptau_pz) private(i) schedule(dynamic, 1)
    for (i = 0; i < N; i++) {
        const VecXd& q_meas = trajPtr_->q(i);
        const VecXd& q_d_meas = trajPtr_->q_d(i);

        // consider sensor noise here
        VecXInt q(q_meas.size());
        VecXInt q_d(q_d_meas.size());
        for (int j = 0; j < q_meas.size(); j++) {
            q(j) = q_meas(j) + IntervalHelper::makeErrorInterval(sensor_noise.position_error(j), 
                                                                 sensor_noise.position_error_type, 
                                                                 q_meas(j));
            q_d(j) = q_d_meas(j) + IntervalHelper::makeErrorInterval(sensor_noise.velocity_error(j), 
                                                                     sensor_noise.velocity_error_type, 
                                                                     q_d_meas(j));
        }

        // below is the extended regressor algorithm from ROAM-Lab
        // refer to https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/dynamics/Hqd_and_CTqd.m

        Vec6Int vJ;
        Mat6Int XJ, dXJdq;
        MatRegressorInt bodyRegressor;
        Vec6Int Sd;
        MatXInt pSd_pz(6, 2 * modelPtr_->nv);

        Eigen::Array<Mat6Int, 1, Eigen::Dynamic> Xup(modelPtr_->nv);
        Eigen::Array<Mat6Int, 1, Eigen::Dynamic> dXupdq(modelPtr_->nv);
        Eigen::Array<Vec6Int, 1, Eigen::Dynamic> S(modelPtr_->nv);
        Eigen::Array<Vec6Int, 1, Eigen::Dynamic> v(modelPtr_->nv);
        Eigen::Array<MatXInt, 1, Eigen::Dynamic> pv_pz(modelPtr_->nv);
        Eigen::Array<MatRegressorInt, 1, Eigen::Dynamic> pbodyRegressor_pz(2 * modelPtr_->nv);

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
                    pv_pz(j) = MatXInt::Zero(6, 2 * modelPtr_->nv);

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
                for (int k = 0; k < 2 * modelPtr_->nv; k++) {
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
                    for (int k = 0; k < 2 * modelPtr_->nv; k++) {
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
                    for (int k = 0; k < 2 * modelPtr_->nv; k++) {
                        pY_CTv_pz(k).block(i * modelPtr_->nv + h, j * 10, 1, 10) = 
                            intervalMatrixMultiply(pSd_pz.col(k).transpose(), bodyRegressor) +
                            intervalMatrixMultiply(Sd.transpose(), pbodyRegressor_pz(k));
                    }
                }

                if (modelPtr_->parents[h + 1] > 0) {
                    if (compute_derivatives) {
                        for (int k = 0; k < 2 * modelPtr_->nv; k++) {
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
            for (int k = 0; k < 2 * modelPtr_->nv; k++) {
                ptau_pz(i).col(k) = intervalDoubleMatrixMultiply(pY_pz(k).middleRows(i * modelPtr_->nv, modelPtr_->nv), phi);
            }
        }
    }

    if (compute_derivatives) {
        mrPtr_->compute(z, false);

        #pragma omp parallel for shared(mrPtr_, tau, ptau_pz, Y, pY_pz) private(i) schedule(dynamic, 1)
        for (i = 0; i < trajPtr_->N; i++) {
            for (int j = 0; j < modelPtr_->nv; j++) {
                const Interval position_error = IntervalHelper::makeErrorInterval(sensor_noise.position_error(j), 
                                                                                  sensor_noise.position_error_type, 
                                                                                  trajPtr_->q(i)(j));
                const Interval velocity_error = IntervalHelper::makeErrorInterval(sensor_noise.velocity_error(j),
                                                                                  sensor_noise.velocity_error_type, 
                                                                                  trajPtr_->q_d(i)(j));

                // an alternative way to compute the bound of tau using first order Taylor expansion
                Interval tau_alt = Interval(mrPtr_->tau(i)(j));
                for (int k = 0; k < modelPtr_->nv; k++) {
                    tau_alt += ptau_pz(i)(j, k) * position_error;
                    tau_alt += ptau_pz(i)(j, k + modelPtr_->nv) * velocity_error;
                }

                // update tau if the alternative tau bound is tighter
                if (IntervalHelper::getRadius(tau_alt) < IntervalHelper::getRadius(tau(i)(j))) {
                    tau(i)(j) = tau_alt;
                }

                for (int h = 0; h < 10 * modelPtr_->nv; h++) {
                    Interval Y_alt = Interval(mrPtr_->Y(i * modelPtr_->nv + j, h));
                    Interval Y_CTv_alt = Interval(mrPtr_->Y_CTv(i * modelPtr_->nv + j, h));

                    for (int k = 0; k < modelPtr_->nv; k++) {
                        Y_alt += pY_pz(k)(i * modelPtr_->nv + j, h) * position_error;
                        Y_CTv_alt += pY_CTv_pz(k)(i * modelPtr_->nv + j, h) * position_error;
                    }
                    for (int k = 0; k < modelPtr_->nv; k++) {
                        Y_alt += pY_pz(k + modelPtr_->nv)(i * modelPtr_->nv + j, h) * velocity_error;
                        Y_CTv_alt += pY_CTv_pz(k + modelPtr_->nv)(i * modelPtr_->nv + j, h) * velocity_error;
                    }

                    if (IntervalHelper::getRadius(Y_alt) < IntervalHelper::getRadius(Y(i * modelPtr_->nv + j, h))) {
                        Y(i * modelPtr_->nv + j, h) = Y_alt;
                    }
                    if (IntervalHelper::getRadius(Y_CTv_alt) < IntervalHelper::getRadius(Y_CTv(i * modelPtr_->nv + j, h))) {
                        Y_CTv(i * modelPtr_->nv + j, h) = Y_CTv_alt;
                    }
                }
            }
        }
    }
}

}; // namespace RAPTOR