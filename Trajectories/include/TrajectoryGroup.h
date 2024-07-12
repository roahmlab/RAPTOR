#ifndef TRAJECTORY_GROUP_H
#define TRAJECTORY_GROUP_H

#include <string>
#include <unordered_map>
#include "Trajectories.h"

namespace IDTO {

class TrajectoryGroup : public Trajectories {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    TrajectoryGroup() = default;

    // Destructor
    ~TrajectoryGroup() = default;

    // class methods:
    virtual void add_trajectory(const std::string& name,    
                                std::shared_ptr<Trajectories> trajectory) final override;

        // update original information in Trajectories class
    virtual void gather_trajectories_information(const bool print_info = false) final override;

    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false);

    // class members:
    std::unordered_map<std::string, std::shared_ptr<Trajectories>> trajectories;

        // we will assemble the trajectories in the order they are added in q, q_d, q_dd.
        // this records the start and the end of each trajectory 
        // in the overall time instances of the assembled trajectory.
        // range is start <= ... < end
    std::unordered_map<std::string, std::pair<size_t, size_t>> trajectory_locations;

        // we will assemble the trajectories in the order they are added in q, q_d, q_dd.
        // this records the start and the end of the decision variable of each trajectory 
        // in the overall decision variable vector z.
        // range is start <= ... < end
    std::unordered_map<std::string, std::pair<size_t, size_t>> variable_locations;
};

}; // namespace IDTO    

#endif // TRAJECTORIES_H
