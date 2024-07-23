#ifndef CONTACT_CONSTRAINTS_H
#define CONTACT_CONSTRAINTS_H

#include "Constraints.h"
#include "Trajectories.h"
#include "CustomizedInverseDynamics.h"
#include "Utils.h"

namespace RAPTOR {

typedef struct contactSurfaceParams_  {
    double mu = 0.7; // friction coefficient
    double Lx = 0.1; // radius of the contact surface (assumed to be rectangle)
    double Ly = 0.1; // radius of the contact surface (assumed to be rectangle)
    double maxSuctionForce = 100; // maximum suction force

    contactSurfaceParams_() = default;

    contactSurfaceParams_(double mu_input, 
                          double Lx_input, 
                          double Ly_input,
                          double maxSuctionForce_input) :
        mu(mu_input), 
        Lx(Lx_input), 
        Ly(Ly_input),
        maxSuctionForce(maxSuctionForce_input) {}
} contactSurfaceParams;

class ContactConstraints : public Constraints {
public:
    using Force = pinocchio::Data::Force;
    using Vec3 = Eigen::Vector3d;
    using Vec6 = Vector6d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    ContactConstraints() = default;

    // Constructor
    ContactConstraints(std::shared_ptr<CustomizedInverseDynamics>& idPtr_input,
                       const contactSurfaceParams& csp_input);

    // Destructor
    ~ContactConstraints() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() override;

    // class members:
    std::shared_ptr<CustomizedInverseDynamics> idPtr_;

    contactSurfaceParams csp;
};

}; // namespace RAPTOR

#endif // CONTACT_CONSTRAINTS_H
