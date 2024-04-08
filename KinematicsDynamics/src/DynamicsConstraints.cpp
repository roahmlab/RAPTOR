#include "DynamicsConstraints.h"

namespace IDTO {

void DynamicsConstraints::reinitialize() {

}

DynamicsConstraints::DynamicsConstraints(const int numJoints_input, int numDependentJoints_input) :
    numJoints(numJoints_input),
    numDependentJoints(numDependentJoints_input) {
    numIndependentJoints = numJoints - numDependentJoints;

    c = VecX::Zero(numJoints);
    J = MatX::Zero(numDependentJoints, numJoints);
    Jx_partial_dq = MatX::Zero(numDependentJoints, numJoints);
    JTx_partial_dq = MatX::Zero(numJoints, numJoints);
    Jxy_partial_dq = MatX::Zero(numDependentJoints, numJoints);

    pq_dep_pq_indep = MatX::Zero(numDependentJoints, numIndependentJoints);
    pv_dep_pq = MatX::Zero(numDependentJoints, numJoints);
    pv_dep_pv_indep = MatX::Zero(numDependentJoints, numIndependentJoints);
    pa_dep_pq = MatX::Zero(numDependentJoints, numJoints);
    pa_dep_pv = MatX::Zero(numDependentJoints, numJoints);
    pa_dep_pa_indep = MatX::Zero(numDependentJoints, numIndependentJoints);

    J_dep = MatX::Zero(numDependentJoints, numDependentJoints);
    J_indep = MatX::Zero(numDependentJoints, numJoints - numDependentJoints);

    P_dep = MatX::Zero(numDependentJoints, numIndependentJoints);
    Pa_indep = VecX::Zero(numJoints);
    temp1 = MatX::Zero(numDependentJoints, numJoints);
    temp2_1 = VecX::Zero(numJoints);
    temp2 = MatX::Zero(numDependentJoints, numJoints);
    temp3 = MatX::Zero(numDependentJoints, numJoints);
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