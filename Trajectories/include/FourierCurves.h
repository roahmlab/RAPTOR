#ifndef FOURIERCURVES_H
#define FOURIERCURVES_H

#include "Trajectories.h"

namespace IDTO {

class FourierCurves : public Trajectories {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using SpaMatX = Eigen::SparseMatrix<double>;

    FourierCurves() = default;

    FourierCurves(const VecX& tspan_input, int Nact_input, int degree_input);

    FourierCurves(double T_input, int N_input, int Nact_input, TimeDiscretization time_discretization, int degree_input);

    void compute(const VecX& z, bool compute_derivatives = true) override;

    int degree = 0; // degree of the Fourier series

    Eigen::VectorXd F;
    Eigen::VectorXd dF;
    Eigen::VectorXd ddF;

    Eigen::VectorXd F0;
    Eigen::VectorXd dF0;
};

}; // namespace IDTO

#endif // FOURIERCURVES_H
