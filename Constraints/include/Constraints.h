
#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <Eigen/Dense>

namespace IDTO {

// This is the base (abstract) class for all constraints
class Constraints {
public:
    using VecX = Eigen::VectorXd;

    // Constructor
    Constraints() = default;

    // Destructor
    ~Constraints() = default;

    // class methods:
    virtual void compute(const VecX& z, bool compute_derivatives = true) = 0;

    virtual int return_m() {return m;}

    virtual void compute_lb() = 0;

    virtual void compute_ub() = 0;

    // class members:
    int m = 0; // number of constraints

    // define the variables that stores the results here
};

}; // namespace IDTO

#endif // CONSTRAINTS_H
