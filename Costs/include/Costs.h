#ifndef COSTS_H
#define COSTS_H

#include <Eigen/Dense>
#include <iostream>
#include "Utils.h"

namespace RAPTOR {

// This is the base (abstract) class for all constraints
class Costs {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    Costs() = default;

    Costs(int varLength);

    // Destructor
    ~Costs() = default;

    // class methods:
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) = 0;

        // determine if the cost is computed before and save the current decision variable
    bool is_computed(const VecX& z, 
                     bool compute_derivatives,
                     bool compute_hessian);

    void initialize_memory(const int varLength);

    // class members:
        // the decision variable that was evaluated last time
    VecX current_z;
    bool if_compute_derivatives = false;
    bool if_compute_hessian = false;

        // compute results are stored here
    double f;
    VecX grad_f;
    MatX hess_f;
};

}; // namespace RAPTOR

#endif // COSTS_H
