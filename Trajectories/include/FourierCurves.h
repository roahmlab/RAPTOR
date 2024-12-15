#ifndef FOURIERCURVES_H
#define FOURIERCURVES_H

#include "Trajectories.h"

namespace RAPTOR {

class FourierCurves : public Trajectories {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    FourierCurves() = default;

    FourierCurves(const VecX& tspan_input, 
                  int Nact_input, 
                  int degree_input,
                  VecX q0_input = VecX::Zero(0),
                  VecX q_d0_input = VecX::Zero(0));

    FourierCurves(double T_input, 
                  int N_input, 
                  int Nact_input, 
                  TimeDiscretization time_discretization, 
                  int degree_input,
                  VecX q0_input = VecX::Zero(0),
                  VecX q_d0_input = VecX::Zero(0));

    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) final override;

    int degree = 0; // degree of the Fourier series

    VecX F;
    VecX dF;
    VecX ddF;

    VecX F0;
    VecX dF0;

    VecX pF_pw;
    VecX pdF_pw;
    VecX pddF_pw;

    VecX pF0_pw;
    VecX pdF0_pw;

    VecX q0;
    VecX q_d0;
    bool optimize_initial_position = true;
    bool optimize_initial_velocity = true;
};

}; // namespace RAPTOR

#endif // FOURIERCURVES_H
