#ifndef MINIMIZE_JERK_H
#define MINIMIZE_JERK_H

#include <memory>

#include "Costs.h"
#include "Trajectories.h"
#include "Utils.h"

namespace RAPTOR {
/*
minimize the torque required over the entire trajectory
*/
class MinimizeJerk : public Costs {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    MinimizeJerk() = default;

    // Constructor
    MinimizeJerk(std::shared_ptr<Trajectories>& trajPtr_input);

    // Destructor
    ~MinimizeJerk() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

    // class members:
    std::shared_ptr<Trajectories> trajPtr_;
};

}; // namespace RAPTOR

#endif // MINIMIZE_JERK_H
