#ifndef MINIMIZE_PATH_LENGTH_H
#define MINIMIZE_PATH_LENGTH_H

#include <memory>

#include "Costs.h"
#include "Trajectories.h"
#include "Utils.h"

namespace RAPTOR {
/*
minimize the torque required over the entire trajectory
*/
class MinimizePathLength : public Costs {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    MinimizePathLength() = default;

    // Constructor
    MinimizePathLength(std::shared_ptr<Trajectories>& trajPtr_input);

    // Destructor
    ~MinimizePathLength() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

    // class members:
    std::shared_ptr<Trajectories> trajPtr_;
};

}; // namespace RAPTOR

#endif // MINIMIZE_PATH_LENGTH_H
