#ifndef CONTACT_CONSTRAINTS_H
#define CONTACT_CONSTRAINTS_H

#include "Constraints.h"
#include "Trajectories.h"
#include "ConstrainedInverseDynamics.h"
#include "Utils.h"

namespace RAPTOR {

typedef struct rectangleContactSurfaceParams_  {
    float mu = 0.7;
    float gamma = 0.7;
    float Lx = 0.1;
    float Ly = 0.1;

    rectangleContactSurfaceParams_(
        float mu_input, 
        float gamma_input, 
        float Lx_input, 
        float Ly_input) :
        mu(mu_input), 
        gamma(gamma_input), 
        Lx(Lx_input), 
        Ly(Ly_input) {}
} rectangleContactSurfaceParams;

class RectangleSurfaceContactConstraints : public Constraints {
public:
    // Constructor
    RectangleSurfaceContactConstraints() = default;

    // Constructor
    RectangleSurfaceContactConstraints(std::shared_ptr<ConstrainedInverseDynamics>& cidPtr_input,
                                       const rectangleContactSurfaceParams& fp_input);

    // Destructor
    ~RectangleSurfaceContactConstraints() = default;

    // class methods:
        // compute constraints
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) final override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() final override;

        // print violation info
    void print_violation_info() final override;

    // class variables:
    std::shared_ptr<ConstrainedInverseDynamics> cidPtr_;

    rectangleContactSurfaceParams fp;
};

}; // namespace RAPTOR

#endif // CONTACT_CONSTRAINTS_H
