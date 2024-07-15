#ifndef CONTACT_CONSTRAINTS_H
#define CONTACT_CONSTRAINTS_H

#include "Constraints.h"
#include "Trajectories.h"
#include "ConstrainedInverseDynamics.h"
#include "Utils.h"

namespace IDTO {

typedef struct frictionParams_  {
    double mu = 0.7;
    double gamma = 0.7;
    double Lx = 0;
    double Ly = 0;

    frictionParams_(double mu_input, 
                    double gamma_input, 
                    double Lx_input, 
                    double Ly_input) :
        mu(mu_input), 
        gamma(gamma_input), 
        Lx(Lx_input), 
        Ly(Ly_input) {}
} frictionParams;

class SurfaceContactConstraints : public Constraints {
public:
    // Constructor
    SurfaceContactConstraints() = default;

    // Constructor
    SurfaceContactConstraints(std::shared_ptr<ConstrainedInverseDynamics>& cidPtr_input,
                              const frictionParams& fp_input);

    // Destructor
    ~SurfaceContactConstraints() = default;

    // class methods:
        // compute constraints
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() override;

    // class variables:
    std::shared_ptr<ConstrainedInverseDynamics> cidPtr_;

    frictionParams fp;
};

}; // namespace IDTO

#endif // CONTACT_CONSTRAINTS_H
