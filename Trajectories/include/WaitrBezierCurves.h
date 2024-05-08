#ifndef WAITR_BEZIER_CURVES_H
#define WAITR_BEZIER_CURVES_H

#include "BezierCurves.h"

namespace IDTO {

/*
This class implements a special type of Bezier curves used in [armour](https://arxiv.org/abs/2301.13308)
But we add one more dimension to model the fixed joint between the end effector and the object being manipulated.
The position, velocity and acceleration of the fixed joint are fixed to 0 all the time.
The reason we do that is pinocchio automatically remove the fixed joints and combine two bodies.
We have to add this extra joint to access the derivatives of the force applied on the fixed joint using pinocchio's functions.
*/

const int WAITR_BEZIER_CURVE_DEGREE = 5;

typedef struct WaitrTrajectoryParameters_ {
    Eigen::VectorXd q0;
    Eigen::VectorXd q_d0;
    Eigen::VectorXd q_dd0;
} WaitrTrajectoryParameters;

class WaitrBezierCurves : public BezierCurves {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    WaitrBezierCurves() = default;

    WaitrBezierCurves(const VecX& tspan_input, 
                      int Nact_input, 
                      const WaitrTrajectoryParameters& atp_input);

    WaitrBezierCurves(double T_input, 
                      int N_input, 
                      int Nact_input, 
                      TimeDiscretization time_discretization, 
                      const WaitrTrajectoryParameters& atp_input); 

    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false);

    WaitrTrajectoryParameters atp;

    MatX coefficients;
};

}; // namespace IDTO

#endif // WAITR_BEZIER_CURVES_H
