#ifndef TORQUE_LIMITS_H
#define TORQUE_LIMITS_H

#include "Constraints.h"
#include "InverseDynamics.h"

namespace IDTO {

class TorqueLimits : public Constraints {
public:
    // Constructor
    TorqueLimits() = default;

    // Constructor
    TorqueLimits(std::shared_ptr<Trajectories>& trajPtr_input, 
                 std::shared_ptr<InverseDynamics> idPtr_input,
                 const VecX& lowerLimits_input, 
                 const VecX& upperLimits_input);

    // Destructor
    ~TorqueLimits() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() override;

        // print violation information
    virtual void print_violation_info() override;

    // class members:
    std::shared_ptr<Trajectories> trajPtr_;
    std::shared_ptr<InverseDynamics> idPtr_;
    
    VecX lowerLimits;
    VecX upperLimits;
};

} // namespace IDTO

#endif // TORQUE_LIMITS_H
