
#ifndef DIGIT_CONSTRAINED_INVERSE_DYNAMICS_H
#define DIGIT_CONSTRAINED_INVERSE_DYNAMICS_H

#include "ConstrainedInverseDynamics.h"

#include "DigitDynamicsConstraints.h"
#include "DigitMultipleStepDynamicsConstraints.h"

namespace IDTO {
namespace Digit {

class DigitConstrainedInverseDynamics : public ConstrainedInverseDynamics {
public:
    using Model = pinocchio::Model;

    // Constructor
    DigitConstrainedInverseDynamics() = default;

    // Constructor (for single step optimization)
    DigitConstrainedInverseDynamics(const Model& model_input, 
                                    std::shared_ptr<Trajectories>& trajPtr_input,
                                    int numDependentJoints_input,
                                    const Eigen::VectorXi& jtype_input, 
                                    char stanceLeg, 
                                    const Transform& stance_foot_T_des_input);

    // Constructor (for multiple step optimization)
    DigitConstrainedInverseDynamics(const Model& model_input, 
                                    std::shared_ptr<Trajectories>& trajPtr_input,
                                    int numDependentJoints_input,
                                    const Eigen::VectorXi& jtype_input, 
                                    char stanceLeg, 
                                    const Transform& stance_foot_T_des_input,
                                    const int N_input,
                                    const int num_steps_input);

    // Destructor
    ~DigitConstrainedInverseDynamics() = default;
};

}; // namespace Digit
}; // namespace IDTO

#endif // DIGIT_CONSTRAINED_INVERSE_DYNAMICS_H
