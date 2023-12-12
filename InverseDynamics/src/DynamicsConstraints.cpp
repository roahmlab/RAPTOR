#include "DynamicsConstraints.h"

namespace IDTO {

void DynamicsConstraints::reinitialize() {

}

DynamicsConstraints::DynamicsConstraints(const Model& model_input, int numDependentJoints_input) :
    numDependentJoints(numDependentJoints_input) {
    modelPtr_ = std::make_unique<Model>(model_input);
    numIndependentJoints = modelPtr_->nv - numDependentJoints;

    c = VecX::Zero(modelPtr_->nv);
    J = MatX::Zero(numDependentJoints, modelPtr_->nv);
    Jx_partial_dq = MatX::Zero(numDependentJoints, modelPtr_->nv);
    JTx_partial_dq = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    Jxy_partial_dq = MatX::Zero(numDependentJoints, modelPtr_->nv);

    pq_dep_pq_indep = MatX::Zero(numDependentJoints, numIndependentJoints);
    pv_dep_pq = MatX::Zero(numDependentJoints, modelPtr_->nv);
    pv_dep_pv_indep = MatX::Zero(numDependentJoints, numIndependentJoints);
    pa_dep_pq = MatX::Zero(numDependentJoints, modelPtr_->nv);
    pa_dep_pv = MatX::Zero(numDependentJoints, modelPtr_->nv);
    pa_dep_pa_indep = MatX::Zero(numDependentJoints, numIndependentJoints);

    J_dep = MatX::Zero(numDependentJoints, numDependentJoints);
    J_indep = MatX::Zero(numDependentJoints, modelPtr_->nv - numDependentJoints);

    P_dep = MatX::Zero(numDependentJoints, numIndependentJoints);
    Pa_indep = VecX::Zero(modelPtr_->nv);
    temp1 = MatX::Zero(numDependentJoints, modelPtr_->nv);
    temp2_1 = VecX::Zero(modelPtr_->nv);
    temp2 = MatX::Zero(numDependentJoints, modelPtr_->nv);
    temp3 = MatX::Zero(numDependentJoints, modelPtr_->nv);
}

void DynamicsConstraints::setupJointPositionVelocityAcceleration(VecX& q, VecX& v, VecX& a, bool compute_derivatives) {
    setupJointPosition(q);

    get_c(q); // c is updated
    get_J(q); // J is updated

    get_dependent_columns(J_dep, J);
    get_independent_columns(J_indep, J);

    J_dep_qr = QRSolver(J_dep);
    J_dep_T_qr = QRSolver(J_dep.transpose());

    // sanity check on uniqueness (these two arguments are actually equivalent)
    if (J_dep_qr.rank()   != J_dep.rows() || 
        J_dep_qr.rank()   != J_dep.cols() ||
        J_dep_T_qr.rank() != J_dep.rows() || 
        J_dep_T_qr.rank() != J_dep.cols()) {
        throw std::runtime_error("constraint jacobian is not full rank!");
    }

    P_dep = -J_dep_qr.solve(J_indep);

    // fill in depuated joints velocities
    fill_dependent_vector(v, P_dep * get_independent_vector(v));

    get_Jx_partial_dq(q, v);
    MatX Jv_partial_dq = Jx_partial_dq;

    // fill in depuated joints accelerations
    fill_dependent_vector(a, -J_dep_qr.solve(J_indep * get_independent_vector(a) + 
                                             Jv_partial_dq * v));              

    if (compute_derivatives) {
        pq_dep_pq_indep = P_dep;

        pv_dep_pq = -J_dep_qr.solve(Jx_partial_dq);
        pv_dep_pv_indep = P_dep;
            
        fill_dependent_vector(Pa_indep, P_dep * get_independent_vector(a));
        fill_independent_vector(Pa_indep, get_independent_vector(a));

        get_Jx_partial_dq(q, Pa_indep); // be careful here, Jx_partial_dq has been changed!
        temp1 = Jx_partial_dq;
        
        fill_dependent_vector(temp2_1, J_dep_qr.solve(Jv_partial_dq * v), true);
        get_Jx_partial_dq(q, temp2_1); // be careful here, Jx_partial_dq has been changed!
        temp2 = Jx_partial_dq;

        get_Jxy_partial_dq(q, v, v);
        temp3 = Jxy_partial_dq;

        pa_dep_pq = J_dep_qr.solve(-temp1 + temp2 - temp3);
        pa_dep_pv = -J_dep_qr.solve(2 * Jv_partial_dq);
        pa_dep_pa_indep = P_dep;
    }      
}

}; // namespace IDTO