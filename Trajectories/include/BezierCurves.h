#ifndef BEZIER_CURVES_H
#define BEZIER_CURVES_H

#include "Trajectories.h"

namespace RAPTOR {

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

    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false);

    int degree = 0; // degree of the Bezier curve

    VecX B;
    VecX dB;
    VecX ddB;

    VecX Bionomials;
};

}; // namespace RAPTOR

#endif // BEZIER_CURVES_H
