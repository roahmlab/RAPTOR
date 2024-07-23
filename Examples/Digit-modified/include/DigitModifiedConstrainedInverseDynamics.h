
#ifndef DIGIT_MODIFIED_CONSTRAINED_INVERSE_DYNAMICS_H
#define DIGIT_MODIFIED_CONSTRAINED_INVERSE_DYNAMICS_H

#include "ConstrainedInverseDynamics.h"
#include "DigitModifiedDynamicsConstraints.h"

namespace RAPTOR {
namespace DigitModified {

class DigitModifiedConstrainedInverseDynamics : public ConstrainedInverseDynamics {
public:
    using Model = pinocchio::Model;

    // Constructor
    DigitModifiedConstrainedInverseDynamics() = default;

    // Constructor
    DigitModifiedConstrainedInverseDynamics(const Model& model_input, 
                                            std::shared_ptr<Trajectories>& trajPtr_input,
                                            int numDependentJoints_input,
                                            const Eigen::VectorXi& jtype_input, 
                                            char stanceLeg, 
                                            const Transform& stance_foot_T_des_input);

    // Destructor
    ~DigitModifiedConstrainedInverseDynamics() = default;

    // class members:
    // a pointer type of DigitModifiedDynamicsConstraints, 
    // that shares the same memory with dcPtr_ defined in base class ConstrainedInverseDynamics
    // so that other Digit-related class can access specific field in DigitModifiedDynamicsConstraints
    // such as stanceLeg, stance_foot_T_des, etc.
    std::shared_ptr<DigitModifiedDynamicsConstraints> ddcPtr_;
};

}; // namespace DigitModified
}; // namespace RAPTOR

#endif // DIGIT_MODIFIED_CONSTRAINED_INVERSE_DYNAMICS_H
