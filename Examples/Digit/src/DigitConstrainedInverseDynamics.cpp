#include "DigitConstrainedInverseDynamics.h"

namespace IDTO {
namespace Digit {

DigitConstrainedInverseDynamics::DigitConstrainedInverseDynamics(const Model& model_input, 
                                                                 int N_input, 
                                                                 int NUM_DEPENDENT_JOINTS_input,
                                                                 std::unique_ptr<DynamicsConstraints>& dynamics_constraints_input) :
    ConstrainedInverseDynamics(model_input, N_input, NUM_DEPENDENT_JOINTS_input, dynamics_constraints_input) {

}

}; // namespace Digit
}; // namespace IDTO