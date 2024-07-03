#include "KinematicsConstraints.h"

namespace IDTO {

namespace LieSpaceResidual {

Eigen::Vector3d translationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr_,
                                    const Eigen::Vector3d& desiredPosition,
                                    Eigen::MatrixXd* gradientPtr_,
                                    Eigen::Array<Eigen::MatrixXd, 3, 1>* hessianPtr_) {
    if (gradientPtr_ != nullptr) {
        *gradientPtr_ = fkPtr_->getTranslationJacobian();
    }

    if (hessianPtr_ != nullptr) {
        fkPtr_->getTranslationHessian(*hessianPtr_);
    }
    
    return fkPtr_->getTranslation() - desiredPosition;
}

Eigen::Vector3d rotationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr_,
                                 const Eigen::Matrix3d& desiredRotation,
                                 Eigen::MatrixXd* gradientPtr_,
                                 Eigen::Array<Eigen::MatrixXd, 3, 1>* hessianPtr_) {
    bool compute_derivatives = (gradientPtr_ != nullptr) ||
                               (hessianPtr_ != nullptr);
    bool compute_hessian = (hessianPtr_ != nullptr);

    // kinematics chain (derivative is only related to these joints)
    const auto& chain = fkPtr_->chain;

    const Eigen::Matrix3d currentRotation = fkPtr_->getRotation();
    Eigen::Matrix3d residualMatrix = desiredRotation.transpose() * currentRotation;

    Eigen::Array<Eigen::Matrix3d, Eigen::Dynamic, 1> dRdq;
    Eigen::Array<Eigen::Matrix3d, Eigen::Dynamic, Eigen::Dynamic> ddRddq;
    if (compute_derivatives) {
        fkPtr_->getRotationJacobian(dRdq);

        for (auto i : chain) {
            dRdq(i) = desiredRotation.transpose() * dRdq(i);
        }

        if (compute_hessian) {
            fkPtr_->getRotationHessian(ddRddq);

            for (auto i : chain) {
                for (auto j : chain) {
                    ddRddq(i, j) = desiredRotation.transpose() * ddRddq(i, j);
                }
            }
        }
    }

    double traceR = residualMatrix.trace();

    Eigen::VectorXd dtraceRdq;
    Eigen::MatrixXd ddtraceRddq;
    if (compute_derivatives) {
        dtraceRdq = Eigen::VectorXd::Zero(dRdq.size());
        for (auto i : chain) {
            dtraceRdq(i) = dRdq(i).trace();
        }

        if (compute_hessian) {
            ddtraceRddq = Eigen::MatrixXd::Zero(ddRddq.rows(), ddRddq.cols());
            for (auto i : chain) {
                for (auto j : chain) {
                    ddtraceRddq(i, j) = ddRddq(i, j).trace();
                }
            }
        }
    }

    double theta = Utils::safeacos((traceR - 1) / 2);

    Eigen::VectorXd dthetadq;
    Eigen::MatrixXd ddthetaddq;
    if (compute_derivatives) {
        dthetadq = Eigen::VectorXd::Zero(dRdq.size());
        for (auto i : chain) {
            dthetadq(i) = 0.5 * safedacosdx((traceR - 1) / 2) * dtraceRdq(i);
        }

        if (compute_hessian) {
            ddthetaddq = Eigen::MatrixXd::Zero(ddRddq.rows(), ddRddq.cols());
            for (auto i : chain) {
                for (auto j : chain) {
                    ddthetaddq(i, j) = 0.5 * 
                        (0.5 * safeddacosddx((traceR - 1) / 2) * dtraceRdq(i) * dtraceRdq(i) +
                         safedacosdx((traceR - 1) / 2)   * ddtraceRddq(i, j));
                }
            }
        }
    }

    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    const Eigen::Matrix3d RRT = residualMatrix - residualMatrix.transpose();
    Eigen::Matrix3d logR = 0.5 * Utils::safexSinx(theta) * RRT;

    if (compute_derivatives) {
        Eigen::MatrixXd& gradient = *gradientPtr_;
        gradient = Eigen::MatrixXd::Zero(3, dRdq.size());

        for (auto i : chain) {
            Eigen::Matrix3d temp1 = 
                Utils::safedxSinxdx(theta) * dthetadq(i) * RRT;
            Eigen::Matrix3d temp2 = 
                Utils::safexSinx(theta) * (dRdq(i) - dRdq(i).transpose());
            gradient.col(i) = Utils::unskew(0.5 * (temp1 + temp2));
        }

        if (compute_hessian) {
            Eigen::Array<Eigen::MatrixXd, 3, 1>& hessian = *hessianPtr_;
            hessian(0) = Eigen::MatrixXd::Zero(ddRddq.rows(), ddRddq.cols());
            hessian(1) = Eigen::MatrixXd::Zero(ddRddq.rows(), ddRddq.cols());
            hessian(2) = Eigen::MatrixXd::Zero(ddRddq.rows(), ddRddq.cols());

            for (auto i : chain) {
                for (auto j : chain) {
                    Eigen::Matrix3d temp1_1 = 
                        Utils::safeddxSinxddx(theta) * dthetadq(i) * dthetadq(j) * RRT;
                    Eigen::Matrix3d temp1_2 = 
                        Utils::safedxSinxdx(theta) * ddthetaddq(i, j) * RRT;
                    Eigen::Matrix3d temp2 = 
                        Utils::safedxSinxdx(theta) * dthetadq(i) * (dRdq(j) - dRdq(j).transpose());
                    Eigen::Matrix3d temp3 = 
                        Utils::safexSinx(theta) * (ddRddq(i, j) - ddRddq(i, j).transpose());
                    Eigen::Vector3d h = Utils::unskew(0.5 * ((temp1_1 + temp1_2) + 2 * temp2 + temp3));

                    hessian(0)(i, j) = h(0);
                    hessian(1)(i, j) = h(1);
                    hessian(2)(i, j) = h(2);
                }
            }
        }
    }

    return Utils::unskew(logR);
}

}; // namespace LieSpaceResidual

KinematicsConstraints::KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                             const Model* model_input,
                                             const Eigen::VectorXi& jtype_input,
                                             const size_t joint_id_input,
                                             const size_t time_id_input,
                                             const Transform& desiredTransform_input,
                                             const Transform endT_input) :
    trajPtr_(trajPtr_input),
    modelPtr_(model_input),
    jtype(jtype_input),
    joint_id(joint_id_input),
    time_id(time_id_input),
    endT(endT_input) {
    fkPtr_ = std::make_unique<ForwardKinematicsSolver>(modelPtr_, jtype);

    if (joint_id > modelPtr_->nq) {
        throw std::invalid_argument("joint_id should not be larger than model.nq");
    }

    if (time_id >= trajPtr_->N) {
        throw std::invalid_argument("time_id should not be larger than number of instances in the trajectory");
    }

    desiredPosition = desiredTransform_input.p;
    desiredRotation = desiredTransform_input.R;

    if (Utils::ifTwoMatrixEqual(desiredRotation, -desiredRotation.transpose())) {
        throw std::invalid_argument("Input matrix is not skew-symmetric");
    }

    constrainPosition = true;
    constrainRotation = true;

    initialize_memory(6, trajPtr_->varLength);
}

KinematicsConstraints::KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                             const Model* model_input,
                                             const Eigen::VectorXi& jtype_input,
                                             const size_t joint_id_input,
                                             const size_t time_id_input,
                                             const Vec3& desiredPosition_input,
                                             const Transform endT_input) :
    trajPtr_(trajPtr_input),
    modelPtr_(model_input),
    jtype(jtype_input),
    joint_id(joint_id_input),
    time_id(time_id_input),
    endT(endT_input) {
    fkPtr_ = std::make_unique<ForwardKinematicsSolver>(modelPtr_, jtype);

    if (joint_id > modelPtr_->nq) {
        throw std::invalid_argument("joint_id should not be larger than model.nq");
    }

    if (time_id >= trajPtr_->N) {
        throw std::invalid_argument("time_id should not be larger than number of instances in the trajectory");
    }

    desiredPosition = desiredPosition_input;

    constrainPosition = true;
    constrainRotation = false;

    initialize_memory(3, trajPtr_->varLength);
}

KinematicsConstraints::KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                             const Model* model_input,
                                             const Eigen::VectorXi& jtype_input,
                                             const size_t joint_id_input,
                                             const size_t time_id_input,
                                             const Mat3& desiredRotation_input,
                                             const Transform endT_input) :
    trajPtr_(trajPtr_input),
    modelPtr_(model_input),
    jtype(jtype_input),
    joint_id(joint_id_input),
    time_id(time_id_input),
    endT(endT_input) {
    fkPtr_ = std::make_unique<ForwardKinematicsSolver>(modelPtr_, jtype);

    if (joint_id > modelPtr_->nq) {
        throw std::invalid_argument("joint_id should not be larger than model.nq");
    }

    if (time_id >= trajPtr_->N) {
        throw std::invalid_argument("time_id should not be larger than number of instances in the trajectory");
    }

    desiredRotation = desiredRotation_input;
    
    if (Utils::ifTwoMatrixEqual(desiredRotation, -desiredRotation.transpose())) {
        throw std::invalid_argument("Input matrix is not skew-symmetric");
    }

    constrainPosition = false;
    constrainRotation = true;

    initialize_memory(3, trajPtr_->varLength);
}

void KinematicsConstraints::compute(const VecX& z, 
                                    bool compute_derivatives,
                                    bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);
    
    const VecX& q = trajPtr_->q(time_id);
    const MatX& pq_pz = trajPtr_->pq_pz(time_id);
    const Eigen::Array<MatX, Eigen::Dynamic, 1>& pq_pz_pz = trajPtr_->pq_pz_pz.col(time_id);

    if (compute_hessian) {
        fkPtr_->compute(0, joint_id, q, &startT, &endT, 2);
    }
    else if (compute_derivatives) {
        fkPtr_->compute(0, joint_id, q, &startT, &endT, 1);
    }
    else {
        fkPtr_->compute(0, joint_id, q, &startT, &endT, 0);
    }

    MatX pg_pq;
    Eigen::Array<MatX, 3, 1> pg_pq_pq;

    if (constrainPosition) {
        if (compute_hessian) {
            g.head(3) = LieSpaceResidual::translationResidual(fkPtr_, desiredPosition, &pg_pq, &pg_pq_pq);
            pg_pz.topRows(3) = pg_pq * pq_pz;
            for (int i = 0; i < 3; i++) {
                pg_pz_pz(i) = pq_pz.transpose() * pg_pq_pq(i) * pq_pz;

                for (int j = 0; j < pq_pz_pz.size(); j++) {
                    pg_pz_pz(i) += pg_pq(i, j) * pq_pz_pz(j);
                }
            }

        }
        else if (compute_derivatives) {
            g.head(3) = LieSpaceResidual::translationResidual(fkPtr_, desiredPosition, &pg_pq);
            pg_pz.topRows(3) = pg_pq * pq_pz;
        }
        else {
            g.head(3) = LieSpaceResidual::translationResidual(fkPtr_, desiredPosition);
        }
    }
    
    if (constrainRotation) {
        if (compute_hessian) {
            g.tail(3) = LieSpaceResidual::rotationResidual(fkPtr_, desiredRotation, &pg_pq, &pg_pq_pq);
            pg_pz.bottomRows(3) = pg_pq * pq_pz;

            int startIdx = constrainPosition ? 3 : 0;

            for (int i = startIdx; i < startIdx + 3; i++) {
                pg_pz_pz(i) = pq_pz.transpose() * pg_pq_pq(i - startIdx) * pq_pz;

                for (int j = 0; j < pq_pz_pz.size(); j++) {
                    pg_pz_pz(i) += pg_pq(i - startIdx, j) * pq_pz_pz(j);
                }
            }
        }
        else if (compute_derivatives) {
            g.tail(3) = LieSpaceResidual::rotationResidual(fkPtr_, desiredRotation, &pg_pq);
            pg_pz.bottomRows(3) = pg_pq * pq_pz;
        }
        else {
            g.tail(3) = LieSpaceResidual::rotationResidual(fkPtr_, desiredRotation);
        }
    }
}

void KinematicsConstraints::compute_bounds() {
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
}

void KinematicsConstraints::print_violation_info() {
    if (constrainPosition) {
        const VecX& g_position = g.head(3);

        if (abs(g_position(0)) > 1e-5) {
            std::cout << "    Error on position x: " 
                      << g_position(0) 
                      << std::endl;
        }
        if (abs(g_position(1)) > 1e-5) {
            std::cout << "    Error on position y: " 
                      << g_position(1) << std::endl;
        }
        if (abs(g_position(2)) > 1e-5) {
            std::cout << "    Error on position z: " 
                      << g_position(2) << std::endl;
        }
    }

    if (constrainRotation) {
        const VecX& g_rotation = g.tail(3);

        if (g_rotation.norm() > 1e-5) {
            std::cout << "    Error on rotation (norm of skew residual): " 
                      << g_rotation.norm() << std::endl;
        }
    }
}

}; // namespace IDTO