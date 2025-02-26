#ifndef MINIMIZE_POWER_H
#define MINIMIZE_POWER_H

#include <memory>

#include "Costs.h"
#include "Trajectories.h"
#include "InverseDynamics.h"
#include "Utils.h"

namespace RAPTOR {

namespace Power {
constexpr double SQUARE_ROOT_THRESHOLD = 1e-8;
}; // namespace Power

/*
minimize the torque required over the entire trajectory
*/
class MinimizePower : public Costs {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    MinimizePower() = default;

    // Constructor
    MinimizePower(std::shared_ptr<Trajectories>& trajPtr_input, 
                  std::shared_ptr<InverseDynamics>& idPtr_input);

    // Destructor
    ~MinimizePower() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

    // class members:
    std::shared_ptr<Trajectories> trajPtr_;
    std::shared_ptr<InverseDynamics> idPtr_;
};

}; // namespace RAPTOR

#endif // MINIMIZE_POWER_H
