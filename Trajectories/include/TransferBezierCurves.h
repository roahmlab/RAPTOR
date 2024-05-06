#ifndef TRANSFER_BEZIER_CURVES_H
#define TRANSFER_BEZIER_CURVES_H

#include "BezierCurves.h"

namespace IDTO {

/*
This class implements a special type of Bezier curves used in [armour](https://arxiv.org/abs/2301.13308)
But we add one more dimension to model the fixed joint between the end effector and the object being manipulated.
The position, velocity and acceleration of the fixed joint are fixed to 0 all the time.
The reason we do that is pinocchio automatically remove the fixed joints and combine two bodies.
We have to add this extra joint to access the derivatives of the force applied on the fixed joint using pinocchio's functions.

Additionally, the initial position and the end position are pre-specified and fixed.
The initial velocity, acceleration and the end velocity, acceleration are fixed to 0.
This fixes 6 degrees of freedom in the trajectory, which are the first 3 and the last 3 Bezier coefficients of the curve.
The rest of the degree of freedom / Bezier coefficients are free to optimize.
*/

typedef struct TransferTrajectoryParameters_ {
    Eigen::VectorXd q0;
    Eigen::VectorXd qT;
} TransferTrajectoryParameters;

class TransferBezierCurves : public BezierCurves {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    TransferBezierCurves() = default;

    TransferBezierCurves(double T_input, 
                         int N_input, 
                         int Nact_input, 
                         const int degree_input,
                         TimeDiscretization time_discretization, 
                         const TransferTrajectoryParameters& ttp_input); 

    void compute(const VecX& z, bool compute_derivatives = true) override;

    TransferTrajectoryParameters ttp;

    MatX coefficients;
};

}; // namespace IDTO

#endif // TRANSFER_BEZIER_CURVES_H
