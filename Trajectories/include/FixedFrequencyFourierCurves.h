#ifndef FIXED_FOURIERCURVES_H
#define FIXED_FOURIERCURVES_H

#include "Trajectories.h"

namespace IDTO {

class FixedFrequencyFourierCurves : public Trajectories {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    FixedFrequencyFourierCurves() = default;

    FixedFrequencyFourierCurves(const VecX& tspan_input, 
                                int Nact_input, 
                                int degree_input,
                                double base_frequency_input = 10);

    FixedFrequencyFourierCurves(double T_input, 
                                int N_input, 
                                int Nact_input, 
                                TimeDiscretization time_discretization, 
                                int degree_input,
                                double base_frequency_input = 10);

    // class methods:
        // compute trajectories
    void compute(const VecX& z, bool compute_derivatives = true) override;

    // class variables:
    double w = 10; // base frequency

    int degree = 0; // degree of the Fourier series

    VecX F;
    VecX dF;
    VecX ddF;

    VecX F0;
    VecX dF0;
};

}; // namespace IDTO

#endif // FIXED_FOURIERCURVES_H