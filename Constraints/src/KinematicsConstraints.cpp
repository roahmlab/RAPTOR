#include "KinematicsConstraints.h"

namespace IDTO {

KinematicsConstraints::KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                               const Model& model_input,
                                               const Eigen::VectorXi& jtype_input,
                                               const size_t joint_id_input,
                                               const size_t time_id_input,
                                               const Transform& desiredTransform_input,
                                               const Transform endT_input) :
    trajPtr_(trajPtr_input),
    jtype(jtype_input),
    joint_id(joint_id_input),
    time_id(time_id_input),
    endT(endT_input) {
    modelPtr_ = std::make_unique<Model>(model_input);
    fkhofPtr_ = std::make_unique<ForwardKinematicsHighOrderDerivative>();

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

    // This is only reserved for position
    jointTJ = MatX::Zero(3, modelPtr_->nv);
    for (int i = 0; i < 3; i++) {
        jointTH(i) = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    }

    initialize_memory(6, trajPtr_->varLength);
}

KinematicsConstraints::KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                               const Model& model_input,
                                               const Eigen::VectorXi& jtype_input,
                                               const size_t joint_id_input,
                                               const size_t time_id_input,
                                               const Vec3& desiredPosition_input,
                                               const Transform endT_input) :
    trajPtr_(trajPtr_input),
    jtype(jtype_input),
    joint_id(joint_id_input),
    time_id(time_id_input),
    endT(endT_input) {
    modelPtr_ = std::make_unique<Model>(model_input);
    fkhofPtr_ = std::make_unique<ForwardKinematicsHighOrderDerivative>();

    if (joint_id > modelPtr_->nq) {
        throw std::invalid_argument("joint_id should not be larger than model.nq");
    }

    if (time_id >= trajPtr_->N) {
        throw std::invalid_argument("time_id should not be larger than number of instances in the trajectory");
    }

    desiredPosition = desiredPosition_input;

    constrainPosition = true;
    constrainRotation = false;

    jointTJ = MatX::Zero(3, modelPtr_->nv);
    for (int i = 0; i < 3; i++) {
        jointTH(i) = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    }

    initialize_memory(3, trajPtr_->varLength);
}

KinematicsConstraints::KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                               const Model& model_input,
                                               const Eigen::VectorXi& jtype_input,
                                               const size_t joint_id_input,
                                               const size_t time_id_input,
                                               const Mat3& desiredRotation_input,
                                               const Transform endT_input) :
    trajPtr_(trajPtr_input),
    jtype(jtype_input),
    joint_id(joint_id_input),
    time_id(time_id_input),
    endT(endT_input) {
    modelPtr_ = std::make_unique<Model>(model_input);
    fkhofPtr_ = std::make_unique<ForwardKinematicsHighOrderDerivative>();

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

    fkhofPtr_->fk(jointT, *modelPtr_, jtype, joint_id, 0, q, endT, startT);

    if (constrainPosition) {
        g.head(3) = jointT.p - desiredPosition;
    }
    
    if (constrainRotation) {
        Mat3 residual = (desiredRotation.transpose() * jointT.R).log();
        g.tail(3) = Utils::unskew(residual);
    }

    if (compute_derivatives) {
        fkhofPtr_->fk_jacobian(dTdq, *modelPtr_, jtype, joint_id, 0, q, endT, startT);
        
        if (constrainPosition) {
            fkhofPtr_->Transform2xyzJacobian(jointTJ, jointT, dTdq);
            pg_pz.topRows(3) = jointTJ * pq_pz;
        }
        
        if (constrainRotation) {
            pg_pz.bottomRows(3).setZero();

            // The following is actually wrong
            // There's no close form solution to the derivative of the log of a matrix
            for (int i = 0; i < trajPtr_->varLength; i++) {
                for (int j = 0; j < dTdq.size(); j++) {
                    Mat3 temp = jointT.R.transpose() * dTdq[j].R * pq_pz(j, i);
                    pg_pz.bottomRows(3).col(i) += Utils::unskew(temp);
                } 
            }
        }
    }

    if (compute_hessian) {
        fkhofPtr_->fk_hessian(ddTddq, *modelPtr_, jtype, joint_id, 0, q, endT, startT);

        if (constrainPosition) {
            fkhofPtr_->Transform2xyzHessian(jointTH, jointT, dTdq, ddTddq);

            // (1) p2_FK_pq2 * pq_pz1 * pq_pz2
            for (int i = 0; i < 3; i++) {
                pg_pz_pz(i) = pq_pz.transpose() * jointTH(i) * pq_pz;
            }
            
            // (2) p_FK_pq * p2_q_pz2
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < trajPtr_->Nact; j++) {
                    pg_pz_pz(i) += jointTJ(i, j) * pq_pz_pz(j);
                }
            }
        }

        if (constrainRotation) {
            // do nothing here
        }
    }
}

void KinematicsConstraints::compute_bounds() {
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
}

}; // namespace IDTO