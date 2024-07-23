#ifndef TRAJECTORIES_H
#define TRAJECTORIES_H

#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <memory>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "Utils.h"

namespace RAPTOR {

enum TimeDiscretization {
    Uniform = 0, 
    Chebyshev
};

class Trajectories {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    Trajectories() = default;

    // Constructor
    Trajectories(const int varLength_input,
                 const VecX& tspan_input, 
                 int Nact_input);

    // Constructor
    Trajectories(const int varLength_input,
                 double T_input, 
                 int N_input, 
                 int Nact_input, 
                 TimeDiscretization time_discretization);

    // Destructor
    ~Trajectories() = default;

    // class methods:
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false);

    bool is_computed(const VecX& z, 
                     bool compute_derivatives,
                     bool compute_hessian);

    // these methods are defined in TrajectoryGroup
    virtual void add_trajectory(const std::string& name,    
                        std::shared_ptr<Trajectories> trajectory){
        throw std::runtime_error("add_trajectory is not implemented in Trajectories class");
    }
    virtual void gather_trajectories_information(const bool print_info = false) {
        throw std::runtime_error("gather_trajectories_information is not implemented in Trajectories class");
    }

    // class members:
    double T = 0; // total time of the trajectory
    int N = 0;    // number of time instances in tspan
    int Nact = 0;   // number of actuated joints in the system
    VecX tspan;   // a column vector of discrete time instances to check constraint violation
    int varLength = 0; // length of the decision variable vector

        // compute results are stored here
    Eigen::Array<VecX, 1, Eigen::Dynamic> q;
    Eigen::Array<VecX, 1, Eigen::Dynamic> q_d;
    Eigen::Array<VecX, 1, Eigen::Dynamic> q_dd;

        // compute results are stored here
    Eigen::Array<MatX, 1, Eigen::Dynamic> pq_pz;
    Eigen::Array<MatX, 1, Eigen::Dynamic> pq_d_pz;
    Eigen::Array<MatX, 1, Eigen::Dynamic> pq_dd_pz;

        // compute results are stored here
    Eigen::Array<MatX, Eigen::Dynamic, Eigen::Dynamic> pq_pz_pz;
    Eigen::Array<MatX, Eigen::Dynamic, Eigen::Dynamic> pq_d_pz_pz;
    Eigen::Array<MatX, Eigen::Dynamic, Eigen::Dynamic> pq_dd_pz_pz;

        // trajectory class is frequently used in the optimization problem
        // so we store the computed results here to avoid recomputation
    VecX current_z;
    bool if_compute_derivatives = false;
    bool if_compute_hessian = false;
};

}; // namespace RAPTOR    

#endif // TRAJECTORIES_H
