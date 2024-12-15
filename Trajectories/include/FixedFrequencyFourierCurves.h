#ifndef FIXED_FOURIERCURVES_H
#define FIXED_FOURIERCURVES_H

#include "Trajectories.h"

namespace RAPTOR {

class FixedFrequencyFourierCurves : public Trajectories {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    FixedFrequencyFourierCurves() = default;

    FixedFrequencyFourierCurves(const VecX& tspan_input, 
                                int Nact_input, 
                                int degree_input,
                                double base_frequency_input = 10,
                                VecX q0_input = VecX::Zero(0),
                                VecX q_d0_input = VecX::Zero(0));

    FixedFrequencyFourierCurves(double T_input, 
                                int N_input, 
                                int Nact_input, 
                                TimeDiscretization time_discretization, 
                                int degree_input,
                                double base_frequency_input = 10,
                                VecX q0_input = VecX::Zero(0),
                                VecX q_d0_input = VecX::Zero(0));

    // class methods:
        // compute trajectories
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) final override;

    // class variables:
    double w = 10; // base frequency

    int degree = 0; // degree of the Fourier series

    VecX F;
    VecX dF;
    VecX ddF;

    VecX F0;
    VecX dF0;

    VecX q0;
    VecX q_d0;
    bool optimize_initial_position = true;
    bool optimize_initial_velocity = true;
};

}; // namespace RAPTOR

#endif // FIXED_FOURIERCURVES_H
