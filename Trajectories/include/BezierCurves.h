#ifndef BEZIER_CURVES_H
#define BEZIER_CURVES_H

#include "Trajectories.h"

namespace IDTO {

class BezierCurves : public Trajectories {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    BezierCurves() = default;

    BezierCurves(const VecX& tspan_input, 
                 int Nact_input, 
                 int degree_input);

    BezierCurves(double T_input, 
                 int N_input, 
                 int Nact_input, 
                 TimeDiscretization time_discretization, 
                 int degree_input);

    void compute(const VecX& z, bool compute_derivatives = true) override;

    int degree = 0; // degree of the Bezier curve

    VecX B;
    VecX dB;
    VecX ddB;

    VecX Bionomials;
};

}; // namespace IDTO

#endif // BEZIER_CURVES_H
