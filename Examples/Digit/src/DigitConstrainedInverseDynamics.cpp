#include "DigitConstrainedInverseDynamics.h"

namespace IDTO {
namespace Digit {

DigitConstrainedInverseDynamics::DigitConstrainedInverseDynamics(const Model& model_input, 
                                                                 int N_input, 
                                                                 int NUM_DEPENDENT_JOINTS_input,
                                                                 const Eigen::VectorXi& jtype_input, 
                                                                 char stanceLeg_input, 
                                                                 const Transform& stance_foot_T_des_input) :
    ConstrainedInverseDynamics(model_input, N_input, NUM_DEPENDENT_JOINTS_input) {
    dynamicsConstraintsPtr_ = std::make_shared<DigitDynamicsConstraints>(model_input, 
                                                                         jtype_input, 
                                                                         stanceLeg_input,
                                                                         stance_foot_T_des_input);
}

}; // namespace Digit
}; // namespace IDTO