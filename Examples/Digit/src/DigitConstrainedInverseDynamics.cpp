#include "DigitConstrainedInverseDynamics.h"

namespace IDTO {
namespace Digit {

DigitConstrainedInverseDynamics::DigitConstrainedInverseDynamics(const Model& model_input, 
                                                                 std::shared_ptr<Trajectories>& trajPtr_input,
                                                                 int NUM_DEPENDENT_JOINTS_input,
                                                                 const Eigen::VectorXi& jtype_input, 
                                                                 char stanceLeg_input, 
                                                                 const Transform& stance_foot_T_des_input) :
    ConstrainedInverseDynamics(model_input, trajPtr_input, NUM_DEPENDENT_JOINTS_input) {
    ddcPtr_ = std::make_shared<DigitDynamicsConstraints>(modelPtr_, 
                                                         jtype_input, 
                                                         stanceLeg_input,
                                                         stance_foot_T_des_input);
    dcPtr_ = ddcPtr_;
}

DigitConstrainedInverseDynamics::DigitConstrainedInverseDynamics(const Model& model_input, 
                                                                 std::shared_ptr<Trajectories>& trajPtr_input,
                                                                 int NUM_DEPENDENT_JOINTS_input,
                                                                 const Eigen::VectorXi& jtype_input, 
                                                                 char stanceLeg_input, 
                                                                 const Transform& stance_foot_T_des_input,
                                                                 const int N_input,
                                                                 const int NSteps_input) :
    ConstrainedInverseDynamics(model_input, trajPtr_input, NUM_DEPENDENT_JOINTS_input) {
    ddcPtr_ = std::make_shared<DigitMultipleStepDynamicsConstraints>(modelPtr_, 
                                                                     jtype_input, 
                                                                     stanceLeg_input,
                                                                     stance_foot_T_des_input,
                                                                     N_input,
                                                                     NSteps_input);
    dcPtr_ = ddcPtr_;
}

}; // namespace Digit
}; // namespace IDTO