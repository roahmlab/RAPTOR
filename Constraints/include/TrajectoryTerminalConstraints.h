#ifndef TRAJECTORY_TERMINAL_CONSTRAINTS_H
#define TRAJECTORY_TERMINAL_CONSTRAINTS_H

#include <memory>

#include "Constraints.h"
#include "Trajectories.h"
#include "Utils.h"

namespace RAPTOR {

class TrajectoryTerminalConstraints : public Constraints {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    TrajectoryTerminalConstraints() = default;

    // Constructor
    TrajectoryTerminalConstraints(std::shared_ptr<Trajectories>& trajPtr_input, 
                                  const VecX desiredPoistion_input = VecX(0),
                                  const VecX desiredVelocity_input = VecX(0),
                                  const VecX desiredAcceleration_input = VecX(0));

    // Destructor
    ~TrajectoryTerminalConstraints() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

        // compute constraints lower bounds and upper bounds
    virtual void compute_bounds() override;

        // print violation information
    virtual void print_violation_info() override;

    // class members:
    std::shared_ptr<Trajectories> trajPtr_;

    VecX desiredPosition;
    VecX desiredVelocity;
    VecX desiredAcceleration;

    bool constrainTerminalPosition = true;
    bool constrainTerminalVelocity = true;
    bool constrainTerminalAcceleration = true;
};

}; // namespace RAPTOR

#endif // TRAJECTORY_TERMINAL_CONSTRAINTS_H
