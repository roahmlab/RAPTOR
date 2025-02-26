#ifndef MINIMIZE_INITIAL_VELOCITY_H
#define MINIMIZE_INITIAL_VELOCITY_H

#include <memory>

#include "Costs.h"
#include "Trajectories.h"
#include "Utils.h"

namespace RAPTOR {

namespace InitialVelocity {
constexpr double SQUARE_ROOT_THRESHOLD = 1e-8;
}; // namespace InitialVelocity

/*
minimize the torque required over the entire trajectory
*/
class MinimizeInitialVelocity : public Costs {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    MinimizeInitialVelocity() = default;

    // Constructor
    MinimizeInitialVelocity(std::shared_ptr<Trajectories>& trajPtr_input);

    // Destructor
    ~MinimizeInitialVelocity() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

    // class members:
    std::shared_ptr<Trajectories> trajPtr_;
};

}; // namespace RAPTOR

#endif // MINIMIZE_INITIAL_VELOCITY_H
