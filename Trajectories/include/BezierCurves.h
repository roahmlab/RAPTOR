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

    void constrainInitialPosition(const VecX& q0_input);

    void constrainInitialVelocity(const VecX& q_d0_input);

    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

    int degree = 0; // degree of the Bezier curve

    VecX B;
    VecX dB;
    VecX ddB;

    VecX Bionomials;

    VecX q0;
    VecX q_d0;
    VecX q_dd0;
    VecX qT;

    bool constrain_initial_position = false;
    bool constrain_initial_velocity = false;
    bool constrain_initial_acceleration = false;
    bool constrain_end_position = false;
};

}; // namespace RAPTOR

#endif // BEZIER_CURVES_H
