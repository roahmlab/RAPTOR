#include "DigitMultipleStepDynamicsConstraints.h"

namespace IDTO {
namespace Digit {

DigitMultipleStepDynamicsConstraints::DigitMultipleStepDynamicsConstraints(const std::shared_ptr<Model>& modelPtr_input,
                                                                           const Eigen::VectorXi& jtype_input, 
                                                                           char stanceLeg, 
                                                                           const Transform& stance_foot_T_des_input,
                                                                           const int N_input,
                                                                           const int NSteps_input) :
    DigitDynamicsConstraints(modelPtr_input, jtype_input, stanceLeg, stance_foot_T_des_input),
    N(N_input),
    NSteps(NSteps_input) {
}

void DigitMultipleStepDynamicsConstraints::setupJointPosition(VecX& q, bool compute_derivatives) {
    DigitDynamicsConstraints::setupJointPosition(q, compute_derivatives);

    // record the number that this function has been called 
    counter1++;

    // reinitialize the counter if it reaches the number of evaluation in one walking step,
    // in other words, it already reaches the next walking step and should switch stance foot
    if (counter1 > 0 && counter1 % N == 0) {
        reinitialize();
    }

    // if (firstCall) {
    //     firstCall = false;
    // }

    // // clear the counter if it reaches the number of evaluation in all walking steps
    // if (counter == N * NSteps) {
    //     counter = 0;
    // }
}

void DigitMultipleStepDynamicsConstraints::setupJointPositionVelocityAcceleration(VecX& q, VecX& v, VecX& a, bool compute_derivatives) {
    if (recoverSavedData(q, v, a, compute_derivatives)) {
        return;
    }

    DigitDynamicsConstraints::setupJointPosition(q, compute_derivatives);

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

    // record the number that this function has been called 
    counter2++;

    // reinitialize the counter if it reaches the number of evaluation in one walking step,
    // in other words, it already reaches the next walking step and should switch stance foot
    if (counter2 > 0 && counter2 % N == 0) {
        reinitialize();
    }

    // if (firstCall) {
    //     firstCall = false;
    // }

    // // clear the counter if it reaches the number of evaluation in all walking steps
    // if (counter == N * NSteps) {
    //     counter = 0;
    // }   
}

}; // namespace Digit
}; // namespace IDTO