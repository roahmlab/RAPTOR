#ifndef MINIMIZE_INITIAL_ACCELERATION_H
#define MINIMIZE_INITIAL_ACCELERATION_H

#include <memory>

#include "Costs.h"
#include "Trajectories.h"
#include "Utils.h"

namespace RAPTOR {

namespace InitialAcceleration {
constexpr double SQUARE_ROOT_THRESHOLD = 1e-8;
}; // namespace InitialAcceleration

/*
minimize the torque required over the entire trajectory
*/
class MinimizeInitialAcceleration : public Costs {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    MinimizeInitialAcceleration() = default;

    // Constructor
    MinimizeInitialAcceleration(std::shared_ptr<Trajectories>& trajPtr_input);

    // Destructor
    ~MinimizeInitialAcceleration() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

    // class members:
    std::shared_ptr<Trajectories> trajPtr_;
};

}; // namespace RAPTOR

#endif // MINIMIZE_INITIAL_ACCELERATION_H
