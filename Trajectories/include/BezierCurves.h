#ifndef BEZIER_CURVES_H
#define BEZIER_CURVES_H

#include "Trajectories.h"

namespace IDTO {

// UN FINISHED

class BezierCurves : public Trajectories {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    BezierCurves() = default;

    BezierCurves(const VecX& tspan_input, int Nact_input, int degree_input);

    BezierCurves(double T_input, int N_input, int Nact_input, TimeDiscretization time_discretization, int degree_input);

    void fixConditionsCheck();

    void fixInitialPosition(const VecX& q0_input);

    void fixInitialVelocity(const VecX& q_d0_input);

    void fixInitialAcceleration(const VecX& q_dd0_input);

    void fixTerminalPosition(const VecX& qf_input);

    void fixTerminalVelocity(const VecX& q_df_input);

    void fixTerminalAcceleration(const VecX& q_ddf_input);

    void compute(const VecX& z, bool compute_derivatives = true) override;

    int degree = 0; // degree of the Bezier curve

    VecX B;
    VecX dB;
    VecX ddB;

    VecX Bionomials;

    // You can use the following variables to fix the initial and terminal conditions
    // so that these conditions are automatically satisfied and there are less decision variables
    VecX q0;
    VecX q_d0;
    VecX q_dd0;
    VecX qf;
    VecX q_df;
    VecX q_ddf;

    bool setInitialPositionFlag = false;
    bool setInitialVelocityFlag = false;
    bool setInitialAccelerationFlag = false;
    bool setTerminalPositionFlag = false;
    bool setTerminalVelocityFlag = false;
    bool setTerminalAccelerationFlag = false;

    int start_idx = 0;
    int end_idx = 0;
};

}; // namespace IDTO

#endif // BEZIER_CURVES_H
