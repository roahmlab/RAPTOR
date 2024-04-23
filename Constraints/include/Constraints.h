
#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <iostream>
#include <Eigen/Dense>
#include <omp.h>
#include "Utils.h"

namespace IDTO {

// This is the base (abstract) class for all constraints
class Constraints {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    Constraints() = default;

    // Destructor
    ~Constraints() = default;

    // class methods:
    virtual void compute(const VecX& z, bool compute_derivatives = true) = 0;

    virtual int return_m() {return m;}

    virtual void compute_bounds() = 0;

    virtual void print_violation_info() {}

        // determine if the constraints are computed before and save the current decision variable
    bool is_computed(const VecX& z, bool compute_derivatives) {
        if (!ifTwoVectorEqual(current_z, z, 0)) {
            current_z = z;
            if_compute_derivatives = compute_derivatives;
            return false;
        }

        if (compute_derivatives != if_compute_derivatives) {
            current_z = z;
            if_compute_derivatives = compute_derivatives;
            return false;
        }

        // current_z = z;  
        if_compute_derivatives = compute_derivatives;
        return true;
    }

    // class members:
    int m = 0; // number of constraints

    // the decision variable that was evaluated last time
    VecX current_z;
    bool if_compute_derivatives = false;

    // compute results are stored here
    VecX g;
    MatX pg_pz;

    VecX g_lb;
    VecX g_ub;

    // define the variables that stores the results here
};

}; // namespace IDTO

#endif // CONSTRAINTS_H
