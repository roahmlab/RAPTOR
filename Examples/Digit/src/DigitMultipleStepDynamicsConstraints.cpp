#include "DigitMultipleStepDynamicsConstraints.h"

namespace IDTO {
namespace Digit {

DigitMultipleStepDynamicsConstraints::DigitMultipleStepDynamicsConstraints(const std::shared_ptr<Model>& modelPtr_input,
                                                                           const Eigen::VectorXi& jtype_input, 
                                                                           char stanceLeg, 
                                                                           const Transform& stance_foot_T_des_input,
                                                                           const int N_input,
                                                                           const int num_steps_input) :
    DigitDynamicsConstraints(modelPtr_input, jtype_input, stanceLeg, stance_foot_T_des_input),
    N(N_input),
    num_steps(num_steps_input) {
}

void DigitMultipleStepDynamicsConstraints::setupJointPosition(VecX& q, bool compute_derivatives) {
    DigitDynamicsConstraints::setupJointPosition(q, compute_derivatives);

    // record the number that this function has been called 
    counter++;

    // reinitialize the counter if it reaches the number of evaluation in one walking step,
    // in other words, it already reaches the next walking step and should switch stance foot
    if (counter == N) {
        reinitialize();
    }

    // clear the counter if it reaches the number of evaluation in all walking steps
    if (counter == N * num_steps) {
        counter = 0;
    }
}

}; // namespace Digit
}; // namespace IDTO