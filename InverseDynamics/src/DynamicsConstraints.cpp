#include "DynamicsConstraints.h"

namespace IDTO {

DynamicsConstraints::DynamicsConstraints(const Model& model_input, int numDependentJoints_input) :
    numDependentJoints(numDependentJoints_input) {
    modelPtr_ = std::make_unique<Model>(model_input);
    numIndependentJoints = modelPtr_->nv - numDependentJoints;

    c = VecX::Zero(modelPtr_->nv);
    J = MatX::Zero(numDependentJoints, modelPtr_->nv);
    Jx_partial_dq = MatX::Zero(numDependentJoints, modelPtr_->nv);
    JTx_partial_dq = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    Jxy_partial_dq = MatX::Zero(numDependentJoints, modelPtr_->nv);

    pq_unact_pq_act = MatX::Zero(numDependentJoints, numIndependentJoints);
    pq_unact_d_pq = MatX::Zero(numDependentJoints, modelPtr_->nv);
    pq_unact_d_pq_act_d = MatX::Zero(numDependentJoints, numIndependentJoints);
    pq_unact_dd_pq = MatX::Zero(numDependentJoints, modelPtr_->nv);
    pq_unact_dd_pq_d = MatX::Zero(numDependentJoints, modelPtr_->nv);
    pq_unact_dd_pq_act_dd = MatX::Zero(numDependentJoints, numIndependentJoints);

    J_dep = MatX::Zero(numDependentJoints, numDependentJoints);
    J_indep = MatX::Zero(numDependentJoints, modelPtr_->nv - numDependentJoints);

    P_dep = MatX::Zero(numDependentJoints, numIndependentJoints);
    Pa_indep = VecX::Zero(modelPtr_->nv);
    temp1 = MatX::Zero(numDependentJoints, modelPtr_->nv);
    temp2_1 = VecX::Zero(modelPtr_->nv);
    temp2 = MatX::Zero(numDependentJoints, modelPtr_->nv);
    temp3 = MatX::Zero(numDependentJoints, modelPtr_->nv);

    v_dep_partial_dq = MatX::Zero(numDependentJoints, modelPtr_->nv);
    a_dep_partial_dq = MatX::Zero(numDependentJoints, modelPtr_->nv);
    a_dep_partial_dv = MatX::Zero(numDependentJoints, modelPtr_->nv);
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

    P_dep = -J_dep_qr.solve(J_indep);

    // fill in unactuated joints velocities
    fill_dependent_vector(v, P_dep * get_independent_vector(v));

    get_Jx_partial_dq(q, v);

    // fill in unactuated joints accelerations
    fill_dependent_vector(a, -J_dep_qr.solve(J_indep * get_independent_vector(a) + 
                                             Jx_partial_dq * v));

    if (compute_derivatives) {
        v_dep_partial_dq = -J_dep_qr.solve(Jx_partial_dq);
            
        fill_dependent_vector(Pa_indep, P_dep * get_independent_vector(a));
        fill_independent_vector(Pa_indep, get_independent_vector(a));

        get_Jx_partial_dq(q, Pa_indep);
        temp1 = Jx_partial_dq;
        
        fill_dependent_vector(temp2_1, J_dep_qr.solve(Jx_partial_dq * v), true);
        get_Jx_partial_dq(q, temp2_1);
        temp2 = Jx_partial_dq;

        get_Jxy_partial_dq(q, v, v);
        temp3 = Jxy_partial_dq;

        a_dep_partial_dq = J_dep_qr.solve(-temp1 + temp2 - temp3);

        a_dep_partial_dv = -J_dep_qr.solve(2 * Jx_partial_dq);

        pq_unact_pq_act = P_dep;
        pq_unact_d_pq = v_dep_partial_dq;
        pq_unact_d_pq_act_d = P_dep;
        pq_unact_dd_pq = a_dep_partial_dq;
        pq_unact_dd_pq_d = a_dep_partial_dv;
        pq_unact_dd_pq_act_dd = P_dep;
    }      
}

}; // namespace IDTO