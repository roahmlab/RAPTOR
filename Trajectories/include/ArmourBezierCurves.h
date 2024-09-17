#ifndef ARMOUR_BEZIER_CURVES_H
#define ARMOUR_BEZIER_CURVES_H

#include "BezierCurves.h"

namespace RAPTOR {

/*
This class implements a special type of Bezier curves used in [armour](https://arxiv.org/abs/2301.13308)
The initial position, velocity and acceleration are fixed and given by the user.    
The final velocity and acceleration are fixed to be 0 (static at the end).
The final position is the trajectory parameter to optimize.
*/

const int ARMOUR_BEZIER_CURVE_DEGREE = 5;

typedef struct ArmourTrajectoryParameters_ {
    Eigen::VectorXf q0;
    Eigen::VectorXf q_d0;
    Eigen::VectorXf q_dd0;
} ArmourTrajectoryParameters;

class ArmourBezierCurves : public BezierCurves {
public:
    using VecX = Eigen::VectorXf;
    using MatX = Eigen::MatrixXf;

    ArmourBezierCurves() = default;

    ArmourBezierCurves(const VecX& tspan_input, 
                       int Nact_input, 
                       const ArmourTrajectoryParameters& atp_input);

    ArmourBezierCurves(float T_input, 
                       int N_input, 
                       int Nact_input, 
                       TimeDiscretization time_discretization, 
                       const ArmourTrajectoryParameters& atp_input); 

    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false);

    ArmourTrajectoryParameters atp;

    MatX coefficients;
};

}; // namespace RAPTOR

#endif // ARMOUR_BEZIER_CURVES_H
