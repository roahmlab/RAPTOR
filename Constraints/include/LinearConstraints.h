#ifndef LINEAR_CONSTRAINTS
#define LINEAR_CONSTRAINTS

#include <memory>

#include "Constraints.h"
#include "Utils.h"

namespace RAPTOR {

class LinearConstraints : public Constraints {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    LinearConstraints() = default;

    // Constructor
    LinearConstraints(const MatX& A_input,
                      const VecX& b_input) :
        Constraints(A_input.rows(), A_input.cols()),
        A(A_input),
        b(b_input) {};

    // Destructor
    ~LinearConstraints() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) {
        if (z.size() != A.cols()) {
            std::cerr << "function input: z.size() = " << z.size() << std::endl;
            std::cerr << "desired: A.cols() = " << A.cols() << std::endl;
            throw std::invalid_argument("LinearConstraints: decision variable vector has wrong size");
        }

        if (is_computed(z, compute_derivatives, compute_hessian)) return;

        g = A * z;

        if (compute_derivatives) {
            pg_pz = A;

            if (compute_hessian) {
                // pg_pz_pz is just zero and has been initialized in the constructor
                // for (int i = 0; i < m; i++) {
                //     pg_pz_pz[i] = MatX::Zero(varLength, varLength);
                // }
            }
        }
    }

        // compute constraints lower bounds and upper bounds
    virtual void compute_bounds() override {
        g_lb = b;
        g_ub = b;
    }

    // class members:
    MatX A;
    VecX b;
};

}; // namespace RAPTOR

#endif // LINEAR_CONSTRAINTS
