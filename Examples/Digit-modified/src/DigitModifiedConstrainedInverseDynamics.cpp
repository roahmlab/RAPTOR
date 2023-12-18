#include "DigitModifiedConstrainedInverseDynamics.h"

namespace IDTO {
namespace DigitModified {

DigitModifiedConstrainedInverseDynamics::DigitModifiedConstrainedInverseDynamics(const Model& model_input, 
                                                                                 std::shared_ptr<Trajectories>& trajPtr_input,
                                                                                 int NUM_DEPENDENT_JOINTS_input,
                                                                                 const Eigen::VectorXi& jtype_input, 
                                                                                 char stanceLeg_input, 
                                                                                 const Transform& stance_foot_T_des_input) :
    ConstrainedInverseDynamics(model_input, trajPtr_input, NUM_DEPENDENT_JOINTS_input) {
    dcPtr_ = std::make_shared<DigitModifiedDynamicsConstraints>(model_input, 
                                                        jtype_input, 
                                                        stanceLeg_input,
                                                        stance_foot_T_des_input);
}

}; // namespace DigitModified
}; // namespace IDTO