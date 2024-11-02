
#ifndef TALOS_CONSTRAINED_INVERSE_DYNAMICS_H
#define TALOS_CONSTRAINED_INVERSE_DYNAMICS_H

#include "ConstrainedInverseDynamics.h"
#include "TalosDynamicsConstraints.h"

namespace RAPTOR {
namespace Talos {

class TalosConstrainedInverseDynamics : public ConstrainedInverseDynamics {
public:
    using Model = pinocchio::Model;

    // Constructor
    TalosConstrainedInverseDynamics() = default;

    // Constructor
    TalosConstrainedInverseDynamics(const Model& model_input, 
                                    std::shared_ptr<Trajectories>& trajPtr_input,
                                    int numDependentJoints_input,
                                    char stanceLeg, 
                                    const Transform& stance_foot_T_des_input);

    // Destructor
    ~TalosConstrainedInverseDynamics() = default;

    // class members:
    // a pointer type of TalosDynamicsConstraints, 
    // that shares the same memory with dcPtr_ defined in base class ConstrainedInverseDynamics
    // so that other Digit-related class can access specific field in TalosDynamicsConstraints
    // such as stanceLeg, stance_foot_T_des, etc.
    std::shared_ptr<TalosDynamicsConstraints> ddcPtr_;
};

}; // namespace Talos
}; // namespace RAPTOR

#endif // TALOS_CONSTRAINED_INVERSE_DYNAMICS_H
