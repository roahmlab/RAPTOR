#include "KinematicsConstraints.h"

namespace IDTO {

KinematicsConstraints::KinematicsConstraints(const Model& model_input,
                                             const Eigen::VectorXi& jtype_input,
                                             std::shared_ptr<Trajectories>& trajPtr_input,
                                             const std::string joint_name_input,
                                             const MatX& lowerLimits_input,
                                             const MatX& upperLimits_input,
                                             const Transform startT_input,
                                             const Transform endT_input,
                                             std::shared_ptr<DynamicsConstraints> dcPtr_input) : 
    trajPtr_(trajPtr_input),
    jtype(jtype_input),
    lowerLimits(lowerLimits_input),
    upperLimits(upperLimits_input),
    startT(startT_input),
    endT(endT_input),
    dcPtr_(dcPtr_input) {
    if (lowerLimits.rows() != 6 || upperLimits.rows() != 6) {
        throw std::invalid_argument("KinematicsConstraints: lowerLimits and upperLimits must have 6 rows (xyz and rpy).");
    }

    if (lowerLimits.cols() != trajPtr_->N || upperLimits.cols() != trajPtr_->N) {
        throw std::invalid_argument("KinematicsConstraints: lowerLimits and upperLimits must have trajPtr_->N columns.");
    }

    modelPtr_ = std::make_unique<Model>(model_input);
    fkhofPtr_ = std::make_unique<ForwardKinematicsHighOrderDerivative>();

    if (modelPtr_->nv > trajPtr_->Nact && dcPtr_ == nullptr) {
        throw std::invalid_argument("KinematicsConstraints: dcPtr_ must be defined if modelPtr_->nv > trajPtr_->Nact, which is a constrained system.");
    }
    else if (modelPtr_->nv == trajPtr_->Nact && dcPtr_ != nullptr) {
        throw std::invalid_argument("KinematicsConstraints: dcPtr_ must be nullptr if modelPtr_->nv == trajPtr_->Nact, which is an unconstrained system.");
    }

    if (modelPtr_->existJointName(joint_name_input)) {
        joint_id = modelPtr_->getJointId(joint_name_input) - 1;
    }
    else {
        throw std::runtime_error("Can not find joint: " + joint_name_input);
    }

    jointTJ = MatX::Zero(6, modelPtr_->nv);
    pq_pz = MatX::Zero(modelPtr_->nv, trajPtr_->varLength);

    m = trajPtr_->N * 6;

    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, trajPtr_->varLength);
}

void KinematicsConstraints::compute(const VecX& z, bool compute_derivatives) {
    if (is_computed(z, compute_derivatives)) {
        return;
    }

    if (compute_derivatives) {
        pg_pz.setZero();
    }

    trajPtr_->compute(z, compute_derivatives);

    // for (int i = 0; i < trajPtr_->N; i++) {
    //     VecX q;
        
    //     if (dcPtr_ == nullptr) { // unconstrained case
    //         q = trajPtr_->q(i);
    //     }
    //     else { // constrained case
    //         q = VecX::Zero(modelPtr_->nq);
    //         dcPtr_->fill_independent_vector(q, trajPtr_->q(i));
    //         dcPtr_->setupJointPosition(q, compute_derivatives);
    //     }

    //     fkhofPtr_->fk(jointT, *modelPtr_, jtype, joint_id, 0, q, endT, startT);

    //     g.block(i * 6, 0, 6, 1) = fkhofPtr_->Transform2xyzrpy(jointT);

    //     if (compute_derivatives) {
    //         fkhofPtr_->fk_jacobian(dTdq, *modelPtr_, jtype, joint_id, 0, q, endT, startT);
    //         fkhofPtr_->Transform2xyzrpyJacobian(jointTJ, jointT, dTdq);

    //         if (dcPtr_ == nullptr) { // unconstrained case
    //             pq_pz = trajPtr_->pq_pz(i);
    //         }
    //         else { // constrained case
    //             // fill in independent joints derivatives directly
    //             for (int j = 0; j < dcPtr_->numIndependentJoints; j++) {
    //                 int indenpendentJointIndex = dcPtr_->return_independent_joint_index(j);
    //                 pq_pz.row(indenpendentJointIndex) = trajPtr_->pq_pz(i).row(j);
    //             }

    //             // compute and fill in dependent joints derivatives
    //             MatX pq_dep_pz = dcPtr_->pq_dep_pq_indep * trajPtr_->pq_pz(i);
    //             for (int j = 0; j < dcPtr_->numDependentJoints; j++) {
    //                 int denpendentJointIndex = dcPtr_->return_dependent_joint_index(j);
    //                 pq_pz.row(denpendentJointIndex) = pq_dep_pz.row(j);
    //             }
    //         }

    //         pg_pz.block(i * 6, 0, 6, trajPtr_->varLength) = jointTJ * pq_pz;
    //     }
    // }
}

void KinematicsConstraints::compute_bounds() {
    for (int i = 0; i < trajPtr_->N; i++) {
        g_lb.block(i * 6, 0, 6, 1) = lowerLimits.col(i);
        g_ub.block(i * 6, 0, 6, 1) = upperLimits.col(i);
    }
}

}; // namespace IDTO