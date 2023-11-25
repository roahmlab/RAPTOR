#include "DynamicsConstraints.h"

namespace IDTO {

DynamicsConstraints::DynamicsConstraints(const Model& model_input, int numDependentJoints) {
    modelPtr_ = std::make_unique<Model>(model_input);

    c = VecX::Zero(modelPtr_->nv);
    J = MatX::Zero(numDependentJoints, modelPtr_->nv);
    Jx_partial_dq = MatX::Zero(numDependentJoints, modelPtr_->nv);
    JTx_partial_dq = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    Jxy_partial_dq = MatX::Zero(numDependentJoints, modelPtr_->nv);

    J_dep = MatX::Zero(numDependentJoints, numDependentJoints);
    J_indep = MatX::Zero(numDependentJoints, modelPtr_->nv - numDependentJoints);
}

void DynamicsConstraints::setupJointPositionVelocityAcceleration(VecX& q, VecX& v, VecX& a, bool compute_derivatives) {
    setupJointPosition(q);

    get_c(q);
    get_J(q);

    get_dependent_columns(J_dep, J);
    get_independent_columns(J_indep, J);

    J_dep_qr = QRSolver(J_dep);
    J_dep_T_qr = QRSolver(J_dep.transpose());

    // sanity check on uniqueness (these two arguments are actually equivalent)
    assert(J_dep_qr.rank() == J_dep.rows() && J_dep_qr.rank() == J_dep.cols());
    assert(J_dep_T_qr.rank() == J_dep.rows() && J_dep_T_qr.rank() == J_dep.cols());

    MatX P_dep = -J_dep_qr.solve(J_indep);

    // fill in unactuated joints velocities
    fill_dependent_vector(v, P_dep * get_independent_vector(v));

    get_Jx_partial_dq(q, v);

    // fill in unactuated joints accelerations
    fill_dependent_vector(a, -J_dep_qr.solve(J_indep * get_independent_vector(a) + 
                                             Jx_partial_dq * v));

    if (compute_derivatives) {

    }      
}

}; // namespace IDTO