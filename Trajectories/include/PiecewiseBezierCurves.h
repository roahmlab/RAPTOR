#ifndef PIECEWISE_BEZIER_CURVES_H
#define PIECEWISE_BEZIER_CURVES_H

#include "BezierCurves.h"

namespace RAPTOR {

class PiecewiseBezierCurves : public BezierCurves {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    PiecewiseBezierCurves() = default;

    PiecewiseBezierCurves(double T_input, 
                          int N_input, 
                          int Nact_input, 
                          int degree_input, 
                          const VecX q0_input = VecX::Zero(0),
                          const VecX qT_input = VecX::Zero(0));

    // Destructor
    ~PiecewiseBezierCurves() = default;

    // class methods:
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) final override;

    // class members:
    MatX coefficients;

    // the initial position of the entire trajectory
    VecX q0;

    // if q0 is not set, 
    // then the initial position is part of the decision variables in the optimization.
    // it will be stored at the end of the decision variable vector z.
    bool optimize_begin_position = false;

    // this indicates the offset of the initial position in the decision variable vector z.
    size_t begin_position_offset = 0;

    // the end position of the entire trajectory
    VecX qT;

    // if qT in TransferTrajectoryParameters is not set,
    // then the end position is part of the decision variables in the optimization.
    // it will be stored at the end of the decision variable vector z.
    bool optimize_end_position = false;

    // this indicates the offset of the end position in the decision variable vector z.
    size_t end_position_offset = 0;
};

}; // namespace RAPTOR

#endif // BEZIER_CURVES_H
