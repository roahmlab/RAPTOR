#include "IntervalRegressorInverseDynamics.h"

namespace RAPTOR {

namespace IntervalHelper {
double getCenter(const Interval& x) {
    return 0.5 * (x.lower() + x.upper());
}
double getRadius(const Interval& x) {
    return 0.5 * (x.upper() - x.lower());
}
Interval makeErrorInterval(const double error, 
                           const SensorNoiseInfo::SensorNoiseType type, 
                           const double value) {
    if (error < 0) {
        throw std::invalid_argument("error should be non-negative!");
    }
    if (type == SensorNoiseInfo::SensorNoiseType::Constant) {
        return Interval(-error, error);
    }
    else { // type == SensorNoiseInfo::SensorNoiseType::Ratio
        return std::abs(value) * Interval(-error, error);
    }
}
}; // namespace IntervalHelper

IntervalRegressorInverseDynamics::IntervalRegressorInverseDynamics(const Model& model_input, 
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
        ptau_pz(i) = MatXInt::Zero(modelPtr_->nv, 3 * modelPtr_->nv);
    }

    Y = MatXInt::Zero(N * modelPtr_->nv, numParams);

    // In the original RegressorInverseDynamics, z represents the decision variable, 
    // such as coefficients of the polynomial trajectory.
    // Here, z represents the q, q_d, q_dd themselves since we are using TrajectoryData class.
    pY_pz.resize(3 * modelPtr_->nv);
    for (int i = 0; i < 3 * modelPtr_->nv; i++) {
        pY_pz(i) = MatXInt::Zero(N * modelPtr_->nv, numParams);
    }

    ridPtr_ = std::make_shared<RegressorInverseDynamics>(*modelPtr_, trajPtr_, false);
}

void IntervalRegressorInverseDynamics::compute(const VecXd& z,
                                               bool compute_derivatives,
                                               bool compute_hessian) {
    if (trajPtr_ == nullptr) {
        throw std::runtime_error("trajPtr_ is not defined yet!");
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    const auto& sensor_noise = trajPtr_->sensor_noise;

    int i = 0;
    #pragma omp parallel for shared(trajPtr_, modelPtr_, jtype, Xtree, a_grav, Y, pY_pz) private(i) schedule(dynamic, 1)
    for (i = 0; i < N; i++) {
        const VecXd& q_meas = trajPtr_->q(i);
        const VecXd& q_d_meas = trajPtr_->q_d(i);
        const VecXd& q_dd_meas = trajPtr_->q_dd(i);

        // consider sensor noise here
        VecXInt q(q_meas.size());
        VecXInt q_d(q_d_meas.size());
        VecXInt q_dd(q_dd_meas.size());
        for (int j = 0; j < q_meas.size(); j++) {
            q(j) = q_meas(j) + IntervalHelper::makeErrorInterval(sensor_noise.position_error(j), 
                                                                 sensor_noise.position_error_type, 
                                                                 q_meas(j));
            q_d(j) = q_d_meas(j) + IntervalHelper::makeErrorInterval(sensor_noise.velocity_error(j), 
                                                                     sensor_noise.velocity_error_type, 
                                                                     q_d_meas(j));
            q_dd(j) = q_dd_meas(j) + IntervalHelper::makeErrorInterval(sensor_noise.acceleration_error(j), 
                                                                       sensor_noise.acceleration_error_type, 
                                                                       q_dd_meas(j));
        }

        // below is the original Roy Featherstone's inverse dynamics algorithm
        // refer to https://royfeatherstone.org/spatial/v2/index.html#ID

        Vec6Int vJ;
        Mat6Int XJ, dXJdq;
        Eigen::Array<Mat6Int, 1, Eigen::Dynamic> Xup(modelPtr_->nv);
        Eigen::Array<Mat6Int, 1, Eigen::Dynamic> dXupdq(modelPtr_->nv);
        Eigen::Array<Vec6Int, 1, Eigen::Dynamic> S(modelPtr_->nv);
        Eigen::Array<Vec6Int, 1, Eigen::Dynamic> v(modelPtr_->nv);
        Eigen::Array<MatXInt, 1, Eigen::Dynamic> pv_pz(modelPtr_->nv);
        Eigen::Array<Vec6Int, 1, Eigen::Dynamic> a(modelPtr_->nv);
        Eigen::Array<MatXInt, 1, Eigen::Dynamic> pa_pz(modelPtr_->nv);

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
                Mat6Int crm_v_j = Spatial::crm(v(j));
                a(j) = intervalMatrixMultiply(Xup(j), a(parent_id)) + intervalMatrixMultiply(crm_v_j, vJ) + S(j) * q_dd(j);

                if (compute_derivatives) {
                    // compute pv_pz
                    pv_pz(j) = intervalMatrixMultiply(Xup(j), pv_pz(parent_id));

                    if (j < trajPtr_->Nact) {
                        // deal with vJ = S(j) * q_d(j);
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

                    // compute pa_pz
                    pa_pz(j) = intervalMatrixMultiply(Xup(j), pa_pz(parent_id));

                    // deal with crm_v_j * vJ;
                    pa_pz(j).row(0) += vJ(2) * pv_pz(j).row(1) - vJ(1) * pv_pz(j).row(2);
                    pa_pz(j).row(1) += vJ(0) * pv_pz(j).row(2) - vJ(2) * pv_pz(j).row(0);
                    pa_pz(j).row(2) += vJ(1) * pv_pz(j).row(0) - vJ(0) * pv_pz(j).row(1);
                    pa_pz(j).row(3) += vJ(2) * pv_pz(j).row(4) - vJ(1) * pv_pz(j).row(5) - vJ(4) * pv_pz(j).row(2) + vJ(5) * pv_pz(j).row(1);
                    pa_pz(j).row(4) += vJ(0) * pv_pz(j).row(5) - vJ(2) * pv_pz(j).row(3) + vJ(3) * pv_pz(j).row(2) - vJ(5) * pv_pz(j).row(0);
                    pa_pz(j).row(5) += vJ(1) * pv_pz(j).row(3) - vJ(0) * pv_pz(j).row(4) - vJ(3) * pv_pz(j).row(1) + vJ(4) * pv_pz(j).row(0);

                    if (j < trajPtr_->Nact) {
                        Vec6Int crm_v_j_S_j = intervalMatrixMultiply(crm_v_j, S(j));
                        for (int k = 0; k < 6; k++) {
                            pa_pz(j)(k, j + modelPtr_->nv) += crm_v_j_S_j(k);
                        }

                        // deal with S(j) * q_dd(j);
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k).lower() != 0 && S(j)(k).upper() != 0) {
                                pa_pz(j)(k, j + 2 * modelPtr_->nv) += S(j)(k);
                            }
                        }

                        // deal with dXupdq(j) * a(parent_id)
                        Vec6Int dXupdq_a_parent_id = intervalMatrixMultiply(dXupdq(j), a(parent_id));
                        for (int k = 0; k < 6; k++) {
                            pa_pz(j)(k, j) += dXupdq_a_parent_id(k);
                        }
                    }
                }
            }
            else {
                v(j) = vJ;
                a(j) = intervalDouble61VectorMultiply(Xup(j), -a_grav) + S(j) * q_dd(j);

                if (compute_derivatives) {
                    // compute pv_pz
                    pv_pz(j) = MatXInt::Zero(6, 3 * modelPtr_->nv);

                    if (j < trajPtr_->Nact) {
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k).lower() != 0 && S(j)(k).upper() != 0) {
                                pv_pz(j)(k, j + modelPtr_->nv) += S(j)(k);
                            }
                        }
                    }

                    // compute pa_pz
                    pa_pz(j) = MatXInt::Zero(6, 3 * modelPtr_->nv);

                    if (j < trajPtr_->Nact) {
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k).lower() != 0 && S(j)(k).upper() != 0) {
                                pa_pz(j)(k, j + 2 * modelPtr_->nv) += S(j)(k);
                            }
                        }

                        Vec6Int dXupdq_a_grav = intervalDouble61VectorMultiply(dXupdq(j), -a_grav);
                        for (int k = 0; k < 6; k++) {
                            pa_pz(j)(k, j) += dXupdq_a_grav(k);
                        }
                    }
                }
            }
        }

        // backward pass
        MatRegressorInt bodyRegressor;
        Eigen::Array<MatRegressorInt, 1, Eigen::Dynamic> pbodyRegressor_pz(3 * modelPtr_->nv);

        for (int j = modelPtr_->nv - 1; j >= 0; j--) {
            // compute local bodyRegressor
            const Interval& v1 = v(j)(0);
            const Interval& v2 = v(j)(1);
            const Interval& v3 = v(j)(2);
            const Interval& v4 = v(j)(3);
            const Interval& v5 = v(j)(4);
            const Interval& v6 = v(j)(5);

            const Interval& a1 = a(j)(0);
            const Interval& a2 = a(j)(1);
            const Interval& a3 = a(j)(2);
            const Interval& a4 = a(j)(3);
            const Interval& a5 = a(j)(4);
            const Interval& a6 = a(j)(5);

            bodyRegressor << 
                                 0,                  0, a6 + v1*v5 - v2*v4, v1*v6 - a5 - v3*v4,     a1,    a2 - v1*v3, -v2*v3,    a3 + v1*v2, v2*v2 - v3*v3,  v2*v3,
                                 0, v2*v4 - v1*v5 - a6,                  0, a4 + v2*v6 - v3*v5,  v1*v3,    a1 + v2*v3,     a2, v3*v3 - v1*v1,    a3 - v1*v2, -v1*v3,
                                 0, a5 - v1*v6 + v3*v4, v3*v5 - v2*v6 - a4,                  0, -v1*v2, v1*v1 - v2*v2,  v1*v2,    a1 - v2*v3,    a2 + v1*v3,     a3,
                a4 + v2*v6 - v3*v5,     -v2*v2 - v3*v3,         v1*v2 - a3,         a2 + v1*v3,      0,             0,      0,             0,             0,      0,
                a5 - v1*v6 + v3*v4,         a3 + v1*v2,     -v1*v1 - v3*v3,         v2*v3 - a1,      0,             0,      0,             0,             0,      0, 
                a6 + v1*v5 - v2*v4,         v1*v3 - a2,         a1 + v2*v3,     -v1*v1 - v2*v2,      0,             0,      0,             0,             0,      0;
                 
            if (compute_derivatives) {
                for (int k = 0; k < 3 * modelPtr_->nv; k++) {
                    const Interval& pv1 = pv_pz(j)(0, k);
                    const Interval& pv2 = pv_pz(j)(1, k);
                    const Interval& pv3 = pv_pz(j)(2, k);
                    const Interval& pv4 = pv_pz(j)(3, k);
                    const Interval& pv5 = pv_pz(j)(4, k);
                    const Interval& pv6 = pv_pz(j)(5, k);

                    const Interval& pa1 = pa_pz(j)(0, k);
                    const Interval& pa2 = pa_pz(j)(1, k);
                    const Interval& pa3 = pa_pz(j)(2, k);
                    const Interval& pa4 = pa_pz(j)(3, k);
                    const Interval& pa5 = pa_pz(j)(4, k);
                    const Interval& pa6 = pa_pz(j)(5, k);

                    pbodyRegressor_pz(k) <<
                                                              0,                                       0, pa6 + v1*pv5 + pv1*v5 - v2*pv4 - pv2*v4, v1*pv6 + pv1*v6 - pa5 - v3*pv4 - pv3*v4,              pa1,   pa2 - v1*pv3 - pv1*v3, -pv2*v3 - v2*pv3,   pa3 + v1*pv2 + pv1*v2, 2.0*v2*pv2 - 2.0*v3*pv3,  pv2*v3 + v2*pv3,
                                                              0, v2*pv4 + pv2*v4 - v1*pv5 - pv1*v5 - pa6,                                       0, pa4 + v2*pv6 + pv2*v6 - v3*pv5 - pv3*v5,  v1*pv3 + pv1*v3,   pa1 + v2*pv3 + pv2*v3,              pa2, 2.0*v3*pv3 - 2.0*v1*pv1,   pa3 - v1*pv2 - pv1*v2, -v1*pv3 - pv1*v3,
                                                              0, v3*pv4 + pv3*v4 - v1*pv6 - pv1*v6 + pa5, v3*pv5 + pv3*v5 - v2*pv6 - pv2*v6 - pa4,                                       0, -v1*pv2 - pv1*v2, 2.0*v1*pv1 - 2.0*v2*pv2,  v1*pv2 + pv1*v2,   pa1 - v2*pv3 - pv2*v3,   pa2 + v1*pv3 + pv1*v3,              pa3,
                        pa4 + v2*pv6 + pv2*v6 - v3*pv5 - pv3*v5,                -2.0*v2*pv2 - 2.0*v3*pv3,                   v1*pv2 + pv1*v2 - pa3,                   pa2 + v1*pv3 + pv1*v3,                0,                       0,                0,                       0,                       0,                0, 
                        v3*pv4 + pv3*v4 - v1*pv6 - pv1*v6 + pa5,                   pa3 + v1*pv2 + pv1*v2,                -2.0*v1*pv1 - 2.0*v3*pv3,                   v2*pv3 + pv2*v3 - pa1,                0,                       0,                0,                       0,                       0,                0,
                        pa6 + v1*pv5 + pv1*v5 - v2*pv4 - pv2*v4,                   v1*pv3 + pv1*v3 - pa2,                   pa1 + v2*pv3 + pv2*v3,                -2.0*v1*pv1 - 2.0*v2*pv2,                0,                       0,                0,                       0,                       0,                0;
                }
            }

            // fill in regressor and backpropagate
            for (int h = j; h >= 0; h = modelPtr_->parents[h + 1] - 1) {
                Y.block(i * modelPtr_->nv + h, j * 10, 1, 10) = 
                    intervalMatrixMultiply(S(j).transpose(), bodyRegressor);

                if (compute_derivatives) {
                    for (int k = 0; k < 3 * modelPtr_->nv; k++) {
                        pY_pz(k).block(i * modelPtr_->nv + h, j * 10, 1, 10) = 
                            intervalMatrixMultiply(S(j).transpose(), pbodyRegressor_pz(k));
                    }
                }

                if (modelPtr_->parents[h + 1] > 0) {
                    if (compute_derivatives) {
                        for (int k = 0; k < 3 * modelPtr_->nv; k++) {
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
            for (int k = 0; k < 3 * modelPtr_->nv; k++) {
                ptau_pz(i).col(k) = intervalDoubleMatrixMultiply(pY_pz(k).middleRows(i * modelPtr_->nv, modelPtr_->nv), phi);
            }
        }
    }

    if (compute_derivatives) {
        ridPtr_->compute(z, false);

        #pragma omp parallel for shared(ridPtr_, tau, ptau_pz, Y, pY_pz) private(i) schedule(dynamic, 1)
        for (i = 0; i < trajPtr_->N; i++) {
            for (int j = 0; j < modelPtr_->nv; j++) {
                const Interval position_error = IntervalHelper::makeErrorInterval(sensor_noise.position_error(j), 
                                                                                  sensor_noise.position_error_type, 
                                                                                  trajPtr_->q(i)(j));
                const Interval velocity_error = IntervalHelper::makeErrorInterval(sensor_noise.velocity_error(j),
                                                                                  sensor_noise.velocity_error_type, 
                                                                                  trajPtr_->q_d(i)(j));
                const Interval acceleration_error = IntervalHelper::makeErrorInterval(sensor_noise.acceleration_error(j),
                                                                                      sensor_noise.acceleration_error_type, 
                                                                                      trajPtr_->q_dd(i)(j));

                // an alternative way to compute the bound of tau using first order Taylor expansion
                Interval tau_alt = Interval(ridPtr_->tau(i)(j));
                for (int k = 0; k < modelPtr_->nv; k++) {
                    tau_alt += ptau_pz(i)(j, k) * position_error;
                    tau_alt += ptau_pz(i)(j, k + modelPtr_->nv) * velocity_error;
                    tau_alt += ptau_pz(i)(j, k + 2 * modelPtr_->nv) * acceleration_error;
                }

                // update tau if the alternative tau bound is tighter
                if (IntervalHelper::getRadius(tau_alt) < IntervalHelper::getRadius(tau(i)(j))) {
                    tau(i)(j) = tau_alt;
                }
                
                for (int h = 0; h < 10 * modelPtr_->nv; h++) {
                    Interval Y_alt = Interval(ridPtr_->Y(i * modelPtr_->nv + j, h));

                    for (int k = 0; k < modelPtr_->nv; k++) {
                        Y_alt += pY_pz(k)(i * modelPtr_->nv + j, h) * position_error;
                    }
                    for (int k = 0; k < modelPtr_->nv; k++) {
                        Y_alt += pY_pz(k + modelPtr_->nv)(i * modelPtr_->nv + j, h) * velocity_error;
                    }
                    for (int k = 0; k < modelPtr_->nv; k++) {
                        Y_alt += pY_pz(k + 2 * modelPtr_->nv)(i * modelPtr_->nv + j, h) * acceleration_error;
                    }

                    if (IntervalHelper::getRadius(Y_alt) < IntervalHelper::getRadius(Y(i * modelPtr_->nv + j, h))) {
                        Y(i * modelPtr_->nv + j, h) = Y_alt;
                    }
                }
            }
        }
    }
}

Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic> intervalMatrixMultiply(
    const Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic>& A,
    const Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic>& B) {
    const size_t RowsA = A.rows();
    const size_t ColsA = A.cols();
    const size_t RowsB = B.rows();
    const size_t ColsB = B.cols();

    if (ColsA != RowsB) {
        throw std::invalid_argument("Matrix dimensions do not match!");
    }

    Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic> result(RowsA, ColsB);

    for (size_t i = 0; i < RowsA; i++) {
        for (size_t j = 0; j < ColsB; j++) {
            result(i, j) = Interval(0); // Initialize interval to zero
            for (int k = 0; k < ColsA; k++) {
                result(i, j) += A(i, k) * B(k, j); // Interval multiplication
            }
        }
    }

    return result;
}

Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic> intervalDoubleMatrixMultiply(
    const Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic>& A,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& B) {
    const size_t RowsA = A.rows();
    const size_t ColsA = A.cols();
    const size_t RowsB = B.rows();
    const size_t ColsB = B.cols();

    if (ColsA != RowsB) {
        throw std::invalid_argument("Matrix dimensions do not match!");
    }

    Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic> result(RowsA, ColsB);

    for (size_t i = 0; i < RowsA; i++) {
        for (size_t j = 0; j < ColsB; j++) {
            result(i, j) = Interval(0); // Initialize interval to zero
            for (int k = 0; k < ColsA; k++) {
                result(i, j) += A(i, k) * B(k, j); // Interval multiplication
            }
        }
    }

    return result;
}

Eigen::Matrix<Interval, 6, 6> intervalDouble66MatrixMultiply(
    const Eigen::Matrix<Interval, 6, 6>& A,
    const Eigen::Matrix<double, 6, 6>& B) {
    Eigen::Matrix<Interval, 6, 6> result;

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            result(i, j) = Interval(0); // Initialize interval to zero
            for (int k = 0; k < 6; k++) {
                result(i, j) += A(i, k) * B(k, j); // Interval multiplication
            }
        }
    }

    return result;
}

Eigen::Matrix<Interval, 6, 1> intervalDouble61VectorMultiply(
    const Eigen::Matrix<Interval, 6, 6>& A,
    const Eigen::Matrix<double, 6, 1>& B) {
    Eigen::Matrix<Interval, 6, 1> result;

    for (int i = 0; i < 6; i++) {
        result(i) = Interval(0); // Initialize interval to zero
        for (int j = 0; j < 6; j++) {
            result(i) += A(i, j) * B(j); // Interval multiplication
        }
    }

    return result;
}

}; // namespace RAPTOR