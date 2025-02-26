#ifndef SUCTION_CUP_CONTACT_CONSTRAINTS_H
#define SUCTION_CUP_CONTACT_CONSTRAINTS_H

#include "Constraints.h"
#include "Trajectories.h"
#include "CustomizedInverseDynamics.h"
#include "Utils.h"

namespace RAPTOR {

typedef struct circleContactSurfaceParams_  {
    double mu = 0.7; // friction coefficient
    double R = 0.1; // radius of the contact surface (assumed to be circle)
    double maxSuctionForce = 0; // maximum suction force of the suction cup

    double contactForceBuffer = 0; // buffer for positive contact force constraint
    double frictionForceBuffer = 0; // buffer for translation friction cone constraint
    double ZMPBuffer = 0; // buffer for ZMP constraint
} circleContactSurfaceParams;

class CircleSurfaceContactConstraints : public Constraints {
public:
    using Force = pinocchio::Data::Force;
    using Vec3 = Eigen::Vector3d;
    using Vec6 = Eigen::Vector<double, 6>;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    CircleSurfaceContactConstraints() = default;

    // Constructor
    CircleSurfaceContactConstraints(std::shared_ptr<CustomizedInverseDynamics>& idPtr_input,
                                    const circleContactSurfaceParams& csp_input);

    // Destructor
    ~CircleSurfaceContactConstraints() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) final override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() final override;

        // print violation info
    void print_violation_info() final override;

    // class members:
    std::shared_ptr<CustomizedInverseDynamics> idPtr_;

    circleContactSurfaceParams csp;
};

}; // namespace RAPTOR

#endif // SUCTION_CUP_CONTACT_CONSTRAINTS_H
