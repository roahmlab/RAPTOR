#ifndef TRAJECTORY_DATA_H
#define TRAJECTORY_DATA_H

#include "Trajectories.h"

namespace RAPTOR {

class TrajectoryData : public Trajectories {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    TrajectoryData() = default;

    TrajectoryData(const std::string& filename_input);

    TrajectoryData(double T_input,
                   int N_input, 
                   int Nact_input);

    // Destructor
    ~TrajectoryData() = default;

    // class methods:
        // has already loaded the data in the constructor
        // simply do nothing here
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) final override {};
};

}; // namespace RAPTOR

#endif // TRAJECTORY_DATA_H
