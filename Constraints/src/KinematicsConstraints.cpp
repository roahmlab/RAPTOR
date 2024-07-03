#include "KinematicsConstraints.h"

namespace IDTO {

namespace LieSpaceResidual {

Eigen::Vector3d translationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr,
                                    const Eigen::Vector3d& desiredPosition) {
    return fkPtr->getTranslation() - desiredPosition;
}

Eigen::Vector3d rotationResidual(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr,
                                 const Eigen::Matrix3d& desiredRotation) {
    const Eigen::Matrix3d currentRotation = fkPtr->getRotation();
    Eigen::Matrix3d residualMatrix = desiredRotation.transpose() * currentRotation;

    // the following is an alternative way to compute the log of a orthogonal matrix
    // Eigen::Matrix3d logR = residualMatrix.log();
    double traceR = residualMatrix.trace();
    double theta = Utils::safeacos((traceR - 1) / 2);
    Eigen::Matrix3d logR = 0.5 * Utils::safexSinx(theta)
                            * (residualMatrix - residualMatrix.transpose());

    // Eigen::Vector3d result;
    // result.setZero();
    // result(0) = traceR;
    // result(1) = theta;
    // result(2) = 0.5 * Utils::safexSinx(theta);
    // return result;
    return Utils::unskew(logR);
}

Eigen::MatrixXd translationResidualGradient(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr,
                                            const Eigen::Vector3d& desiredPosition) {
    return fkPtr->getTranslationJacobian();
}

Eigen::MatrixXd rotationResidualGradient(const std::unique_ptr<ForwardKinematicsSolver>& fkPtr,
                                         const Eigen::Matrix3d& desiredRotation) {
    const Eigen::Matrix3d currentRotation = fkPtr->getRotation();
    Eigen::Matrix3d residualMatrix = desiredRotation.transpose() * currentRotation;

    Eigen::Array<Eigen::Matrix3d, Eigen::Dynamic, 1> dRdq;
    fkPtr->getRotationJacobian(dRdq);

    for (int i = 0; i < dRdq.size(); i++) {
        dRdq(i) = desiredRotation.transpose() * dRdq(i);
    }

    double traceR = residualMatrix.trace();
    Eigen::VectorXd dtraceRdq = Eigen::VectorXd::Zero(dRdq.size());
    for (int i = 0; i < dRdq.size(); i++) {
        dtraceRdq(i) = dRdq(i).trace();
    }

    double theta = Utils::safeacos((traceR - 1) / 2);
    Eigen::VectorXd dthetadq = Eigen::VectorXd::Zero(dRdq.size());
    for (int i = 0; i < dRdq.size(); i++) {
        dthetadq(i) = 0.5 * Utils::safedacosdx((traceR - 1) / 2) * dtraceRdq(i);
    }

    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    const Eigen::Matrix3d RRT = residualMatrix - residualMatrix.transpose();
    Eigen::Matrix3d logR = 0.5 * Utils::safexSinx(theta) * RRT;

    Eigen::MatrixXd result(3, dRdq.size());
    for (int i = 0; i < dRdq.size(); i++) {
        Eigen::Matrix3d temp1 = 
            0.5 * Utils::safedxSinxdx(theta) * dthetadq(i) * RRT;
        Eigen::Matrix3d temp2 = 
            0.5 * Utils::safexSinx(theta) * (dRdq(i) - dRdq(i).transpose());
        result.col(i) = Utils::unskew(temp1 + temp2);
        // result(0, i) = dtraceRdq(i);
        // result(1, i) = dthetadq(i);
        // result(2, i) = 0.5 * Utils::safedxSinxdx(theta) * dthetadq(i);
    }

    return result;
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

    if (constrainPosition) {
        g.head(3) = LieSpaceResidual::translationResidual(fkPtr_, desiredPosition);
    }
    
    if (constrainRotation) {
        g.tail(3) = LieSpaceResidual::rotationResidual(fkPtr_, desiredRotation);
    }

    if (compute_derivatives) {
        if (constrainPosition) {
            const MatX& pg_pq = LieSpaceResidual::translationResidualGradient(fkPtr_, desiredPosition);
            pg_pz.topRows(3) = pg_pq * pq_pz;
        }
        
        if (constrainRotation) {
            const MatX& pg_pq = LieSpaceResidual::rotationResidualGradient(fkPtr_, desiredRotation);
            pg_pz.bottomRows(3) = pg_pq * pq_pz;
        }
    }

    if (compute_hessian) {
        throw std::invalid_argument("KinematicsConstraints: Hessian computation not implemented");
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
        const VecX& g_rotation = g.tail(9);

        if (g_rotation.norm() > 1e-5) {
            std::cout << "    Error on rotation (norm of rotation matrix): " 
                      << g_rotation.norm() << std::endl;
        }
    }
}

}; // namespace IDTO