
#ifndef G1_CONSTRAINED_INVERSE_DYNAMICS_H
#define G1_CONSTRAINED_INVERSE_DYNAMICS_H

#include "ConstrainedInverseDynamics.h"
#include "G1DynamicsConstraints.h"

namespace RAPTOR {
namespace G1 {

class G1ConstrainedInverseDynamics : public ConstrainedInverseDynamics {
public:
    using Model = pinocchio::Model;

    // Constructor
    G1ConstrainedInverseDynamics() = default;

    // Constructor
    G1ConstrainedInverseDynamics(const Model& model_input, 
                                 std::shared_ptr<Trajectories>& trajPtr_input,
                                 int numDependentJoints_input,
                                 char stanceLeg, 
                                 const Transform& stance_foot_T_des_input);

    // Destructor
    ~G1ConstrainedInverseDynamics() = default;

    // class members:
    // a pointer type of G1DynamicsConstraints, 
    // that shares the same memory with dcPtr_ defined in base class ConstrainedInverseDynamics
    // so that other Digit-related class can access specific field in G1DynamicsConstraints
    // such as stanceLeg, stance_foot_T_des, etc.
    std::shared_ptr<G1DynamicsConstraints> ddcPtr_;
};

}; // namespace G1
}; // namespace RAPTOR

#endif // G1_CONSTRAINED_INVERSE_DYNAMICS_H
