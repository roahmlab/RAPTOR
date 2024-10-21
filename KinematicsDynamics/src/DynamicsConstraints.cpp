#include "DynamicsConstraints.h"

namespace RAPTOR {

DynamicsConstraints::DynamicsConstraints(const int numJoints_input, int numDependentJoints_input) :
    numJoints(numJoints_input),
    numDependentJoints(numDependentJoints_input) {
    numIndependentJoints = numJoints - numDependentJoints;

    c = VecX::Zero(numDependentJoints);
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

    bufferDataQueue.resize(0);
}

void DynamicsConstraints::setupJointPositionVelocityAcceleration(VecX& q, VecX& v, VecX& a, bool compute_derivatives) {
    if (recoverSavedData(q, v, a, compute_derivatives)) {
        return;
    }

    setupJointPosition(q, compute_derivatives);

    if (!compute_derivatives) {
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

        // get the dependent columns of the projection matrix
        P_dep = -J_dep_qr.solve(J_indep);
    }
    // else {
    //     // should have been called in setupJointPosition
    // }

    // fill in dependent joints velocities
    fill_dependent_vector(v, P_dep * get_independent_vector(v));

    get_Jx_partial_dq(q, v);
    MatX Jv_partial_dq = Jx_partial_dq;

    // fill in dependent joints accelerations
    fill_dependent_vector(a, -J_dep_qr.solve(J_indep * get_independent_vector(a) + 
                                             Jv_partial_dq * v));              

    if (compute_derivatives) {
        pq_dep_pq_indep = P_dep;

        pv_dep_pq = -J_dep_qr.solve(Jx_partial_dq);
        pv_dep_pv_indep = P_dep;
            
        fill_dependent_vector(Pa_indep, P_dep * get_independent_vector(a));
        fill_independent_vector(Pa_indep, get_independent_vector(a));

        get_Jx_partial_dq(q, Pa_indep); // be careful here, Jx_partial_dq has been changed here!
        temp1 = Jx_partial_dq;
        
        fill_dependent_vector(temp2_1, J_dep_qr.solve(Jv_partial_dq * v), true);
        get_Jx_partial_dq(q, temp2_1); // be careful here, Jx_partial_dq has been changed here!
        temp2 = Jx_partial_dq;

        get_Jxy_partial_dq(q, v, v);
        temp3 = Jxy_partial_dq;

        pa_dep_pq = J_dep_qr.solve(-temp1 + temp2 - temp3);
        pa_dep_pv = -J_dep_qr.solve(2 * Jv_partial_dq);
        pa_dep_pa_indep = P_dep;
    }

    updateQueue(q, v, a, compute_derivatives);      
}

bool DynamicsConstraints::recoverSavedData(VecX& q, bool compute_derivatives) {
    for (auto& it : bufferDataQueue) {
        if (ifIndependentPartEqual(it.q, q) &&
            (it.compute_derivatives == compute_derivatives || it.compute_derivatives == true)) {
            q = it.q;

            if (compute_derivatives) {
                J_dep = it.J_dep;
                J_indep = it.J_indep;
                J_dep_qr = it.J_dep_qr;
                J_dep_T_qr = it.J_dep_T_qr;

                pq_dep_pq_indep = it.pq_dep_pq_indep;
            }

            return true;
        }
    }

    return false;
}

bool DynamicsConstraints::recoverSavedData(VecX& q, VecX& v, VecX& a, bool compute_derivatives) {
    for (auto& it : bufferDataQueue) {
        if (ifIndependentPartEqual(it.q, q) &&
            ifIndependentPartEqual(it.v, v) &&
            ifIndependentPartEqual(it.a, a) &&
            (it.compute_derivatives == compute_derivatives || it.compute_derivatives == true)) {
            q = it.q;
            v = it.v;
            a = it.a;

            J_dep = it.J_dep;
            J_indep = it.J_indep;
            J_dep_qr = it.J_dep_qr;
            J_dep_T_qr = it.J_dep_T_qr;

            if (compute_derivatives) {
                pq_dep_pq_indep = it.pq_dep_pq_indep;
                pv_dep_pq = it.pv_dep_pq;
                pv_dep_pv_indep = it.pv_dep_pv_indep;
                pa_dep_pq = it.pa_dep_pq;
                pa_dep_pv = it.pa_dep_pv;
                pa_dep_pa_indep = it.pa_dep_pa_indep;
            }

            return true;
        }
    }

    return false;
}

void DynamicsConstraints::updateQueue(const VecX& q, bool compute_derivatives) {
    BufferData newData;
    newData.q = q;
    newData.compute_derivatives = compute_derivatives;

    if (compute_derivatives) {
        newData.J_dep = J_dep;
        newData.J_indep = J_indep;
        newData.J_dep_qr = J_dep_qr;
        newData.J_dep_T_qr = J_dep_T_qr;

        newData.pq_dep_pq_indep = pq_dep_pq_indep;
    }

    if (bufferDataQueue.size() > MAX_BUFFER_SIZE) {
        bufferDataQueue.pop_front();
    }

    bufferDataQueue.push_back(newData);
}

void DynamicsConstraints::updateQueue(const VecX& q, const VecX& v, const VecX& a, bool compute_derivatives) {
    BufferData newData;
    newData.q = q;
    newData.v = v;
    newData.a = a;
    newData.compute_derivatives = compute_derivatives;

    newData.J_dep = J_dep;
    newData.J_indep = J_indep;
    newData.J_dep_qr = J_dep_qr;
    newData.J_dep_T_qr = J_dep_T_qr;

    if (compute_derivatives) {
        newData.pq_dep_pq_indep = pq_dep_pq_indep;
        newData.pv_dep_pq = pv_dep_pq;
        newData.pv_dep_pv_indep = pv_dep_pv_indep;
        newData.pa_dep_pq = pa_dep_pq;
        newData.pa_dep_pv = pa_dep_pv;
        newData.pa_dep_pa_indep = pa_dep_pa_indep;
    }

    if (bufferDataQueue.size() > MAX_BUFFER_SIZE) {
        bufferDataQueue.pop_front();
    }

    bufferDataQueue.push_back(newData);
}

}; // namespace RAPTOR