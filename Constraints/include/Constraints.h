#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <Eigen/Dense>
#include <iostream>
#include "Utils.h"

namespace RAPTOR {

// This is the base (abstract) class for all constraints
class Constraints {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    Constraints() = default;

    Constraints(int m_input,
                int varLength);

    // Destructor
    ~Constraints() = default;

    // class methods:
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) = 0;

    virtual void compute_bounds() = 0;

    virtual void print_violation_info() {}

        // determine if the constraints are computed before and save the current decision variable
    bool is_computed(const VecX& z, 
                     bool compute_derivatives,
                     bool compute_hessian);

    void initialize_memory(const int m_input, 
                           const int varLength,
                           bool initialize_hessian = true);

    // class members:
    int m = 0; // number of constraints

    // the decision variable that was evaluated last time
    VecX current_z;
    bool if_compute_derivatives = false;
    bool if_compute_hessian = false;

    // compute results are stored here
    VecX g;
    MatX pg_pz;
    Eigen::Array<MatX, 1, Eigen::Dynamic> pg_pz_pz;

    VecX g_lb;
    VecX g_ub;
};

}; // namespace RAPTOR

#endif // CONSTRAINTS_H
