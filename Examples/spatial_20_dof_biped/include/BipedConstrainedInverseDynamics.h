
#ifndef Biped_CONSTRAINED_INVERSE_DYNAMICS_H
#define Biped_CONSTRAINED_INVERSE_DYNAMICS_H

#include "ConstrainedInverseDynamics.h"
#include "BipedDynamicsConstraints.h"

namespace RAPTOR {
namespace Biped {

class BipedConstrainedInverseDynamics : public ConstrainedInverseDynamics {
public:
    using Model = pinocchio::Model;

    // Constructor
    BipedConstrainedInverseDynamics() = default;

    // Constructor
    BipedConstrainedInverseDynamics(const Model& model_input, 
                                    std::shared_ptr<Trajectories>& trajPtr_input,
                                    int numDependentJoints_input,
                                    char stanceLeg, 
                                    const Transform& stance_foot_T_des_input);

    // Destructor
    ~BipedConstrainedInverseDynamics() = default;

    // class members:
    // a pointer type of BipedDynamicsConstraints, 
    // that shares the same memory with dcPtr_ defined in base class ConstrainedInverseDynamics
    // so that other Digit-related class can access specific field in BipedDynamicsConstraints
    // such as stanceLeg, stance_foot_T_des, etc.
    std::shared_ptr<BipedDynamicsConstraints> ddcPtr_;
};

}; // namespace Biped
}; // namespace RAPTOR

#endif // Biped_CONSTRAINED_INVERSE_DYNAMICS_H
