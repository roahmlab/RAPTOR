#include "EndEffectorConstraints.h"

namespace IDTO {

EndEffectorConstraints::EndEffectorConstraints(const Model& model_input,
                                               const Eigen::VectorXi& jtype_input,
                                               const Transform& endT_input,
                                               const std::string joint_name_input,
                                               std::shared_ptr<Trajectories>& trajPtr_input,
                                               const VecX& desiredEndEffectorPos_input) :
    
    jtype(jtype_input),
    endT(endT_input),
    trajPtr_(trajPtr_input),
    desiredEndEffectorPos(desiredEndEffectorPos_input) {
    if (desiredEndEffectorPos.rows() != 6) {
        throw std::invalid_argument("KinematicsConstraints: desiredEndEffectorPos must have 6 rows (xyz and rpy).");
    }

    modelPtr_ = std::make_unique<Model>(model_input);
    fkhofPtr_ = std::make_unique<ForwardKinematicsHighOrderDerivative>();

    if (modelPtr_->existJointName(joint_name_input)) {
        joint_id = modelPtr_->getJointId(joint_name_input) - 1;
    }
    else {
        throw std::runtime_error("Can not find joint: " + joint_name_input);
    }

    jointTJ = MatX::Zero(6, modelPtr_->nv);
    pq_pz = MatX::Zero(modelPtr_->nv, trajPtr_->varLength);

    m = 6;

    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, trajPtr_->varLength);
}

void EndEffectorConstraints::compute(const VecX& z, bool compute_derivatives) {
    trajPtr_->compute(z, compute_derivatives);

    const VecX& q = trajPtr_->q(trajPtr_->N - 1);

    fkhofPtr_->fk(jointT, *modelPtr_, jtype, joint_id, 0, q, endT, startT);

    g = fkhofPtr_->Transform2xyzrpy(jointT);

    if (compute_derivatives) {
        fkhofPtr_->fk_jacobian(dTdq, *modelPtr_, jtype, joint_id, 0, q, endT, startT);
        fkhofPtr_->Transform2xyzrpyJacobian(jointTJ, jointT, dTdq);

        pq_pz = trajPtr_->pq_pz(trajPtr_->N - 1);
        pg_pz = jointTJ * pq_pz;
    }
}

void EndEffectorConstraints::compute_bounds() {
    g_lb = desiredEndEffectorPos;
    g_ub = desiredEndEffectorPos;
}

}; // namespace IDTO