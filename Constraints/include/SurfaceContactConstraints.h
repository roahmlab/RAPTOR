#ifndef CONTACT_CONSTRAINTS_H
#define CONTACT_CONSTRAINTS_H

#include "Constraints.h"
#include "Trajectories.h"
#include "ConstrainedInverseDynamics.h"

namespace IDTO {

class SurfaceContactConstraints : public Constraints {
public:
    // Constructor
    SurfaceContactConstraints() = default;

    // Constructor
    SurfaceContactConstraints(std::shared_ptr<ConstrainedInverseDynamics>& idPtr_input,
                              const double friction_mu_input,
                              const double friction_gamma_input, 
                              const double Lx_input,
                              const double Ly_input);

    // Destructor
    ~SurfaceContactConstraints() = default;

    // class methods:
        // compute constraints
    void compute(const VecX& z, bool compute_derivatives = true) override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() override;

    // class variables:
    std::shared_ptr<ConstrainedInverseDynamics> idPtr_;

    double friction_mu = 0;
    double friction_gamma = 0;
    double Lx = 0;
    double Ly = 0;
    
};

}; // namespace IDTO

#endif // CONTACT_CONSTRAINTS_H
