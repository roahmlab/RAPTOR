#ifndef TRAJECTORIES_H
#define TRAJECTORIES_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <vector>

namespace IDTO {

enum TimeDiscretization {
    Uniform, 
    Chebyshev
};

class Trajectories {
public:
    using VecX = Eigen::VectorXd;
    using SpaMatX = Eigen::SparseMatrix<double>;

    // Constructor
    Trajectories() = default;

    // Constructor
    Trajectories(const VecX& tspan_input, int Nact_input);

    // Constructor
    Trajectories(double T_input, int N_input, int Nact_input, TimeDiscretization time_discretization);

    // Destructor
    ~Trajectories() = default;

    // class methods:
    virtual void compute(const VecX& z, bool compute_derivatives = true);

    // class members:
    double T = 0; // total time of the trajectory
    int N = 0;    // number of time instances in tspan
    int Nact = 0;   // number of actuated joints in the system
    VecX tspan;   // a column vector of discrete time instances to check constraint violation

        // compute results are stored here
    Eigen::Array<VecX, 1, Eigen::Dynamic> q;
    Eigen::Array<VecX, 1, Eigen::Dynamic> q_d;
    Eigen::Array<VecX, 1, Eigen::Dynamic> q_dd;

        // compute results are stored here
    Eigen::Array<SpaMatX, 1, Eigen::Dynamic> pq_pz;
    Eigen::Array<SpaMatX, 1, Eigen::Dynamic> pq_d_pz;
    Eigen::Array<SpaMatX, 1, Eigen::Dynamic> pq_dd_pz;

        // temporary variables updated in compute()
    std::vector<Eigen::Triplet<double>> pq_pz_tripletList;
    std::vector<Eigen::Triplet<double>> pq_d_pz_tripletList;
    std::vector<Eigen::Triplet<double>> pq_dd_pz_tripletList;
};

}; // namespace IDTO    

#endif // TRAJECTORIES_H
