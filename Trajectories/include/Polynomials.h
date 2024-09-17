#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

#include "Trajectories.h"

namespace RAPTOR {

class Polynomials : public Trajectories {
public:
    using VecX = Eigen::VectorXf;
    using MatX = Eigen::MatrixXf;

    Polynomials() = default;

    Polynomials(const VecX& tspan_input, 
                int Nact_input, 
                int degree_input);

    Polynomials(float T_input,
                int N_input, 
                int Nact_input, 
                TimeDiscretization time_discretization, 
                int degree_input);

    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false);

    int degree = 0; // degree of the polynomial

    VecX P;
    VecX dP;
    VecX ddP;
};

}; // namespace RAPTOR

#endif // POLYNOMIALS_H
