#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

#include "Trajectories.h"

namespace IDTO {

class Polynomials : public Trajectories {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    Polynomials() = default;

    Polynomials(const VecX& tspan_input, 
                int Nact_input, 
                int degree_input);

    Polynomials(double T_input, 
                int N_input, 
                int Nact_input, 
                TimeDiscretization time_discretization, 
                int degree_input);

    void compute(const VecX& z, bool compute_derivatives = true) override;

    int degree = 0; // degree of the polynomial

    VecX P;
    VecX dP;
    VecX ddP;
};

}; // namespace IDTO

#endif // POLYNOMIALS_H
