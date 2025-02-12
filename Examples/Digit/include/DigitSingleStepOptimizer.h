#ifndef DIGIT_SINGLE_STEP_OPTIMIZER_H
#define DIGIT_SINGLE_STEP_OPTIMIZER_H

#include "Optimizer.h"

#include "BezierCurves.h"

#include "DigitConstrainedInverseDynamics.h"
#include "DigitDynamicsConstraints.h"
#include "Utils.h"

#include "ConstrainedJointLimits.h"
#include "TorqueLimits.h"
#include "RectangleSurfaceContactConstraints.h"
#include "DigitCustomizedConstraints.h"
#include "DigitSingleStepPeriodicityConstraints.h"

// #include "MinimizeTorque.h"
#include "MinimizePower.h"
#include "MinimizeInitialVelocity.h"
#include "MinimizeInitialAcceleration.h"

namespace RAPTOR {
namespace Digit {

using namespace Ipopt;

class DigitSingleStepOptimizer : public Optimizer {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    DigitSingleStepOptimizer() = default;

    /** Default destructor */
    ~DigitSingleStepOptimizer() = default;

    // [set_parameters]
    bool set_parameters(
        const VecX& x0_input,
        const double T_input,
        const int N_input,
        const TimeDiscretization time_discretization_input,
        const int degree_input,
        const Model& model_input, 
        const GaitParameters& gp_input,
        const char stanceLeg = 'L', // stance foot is left foot by default
        const Transform& stance_foot_T_des = Transform(3, -M_PI_2),
        bool periodic = true,
        const VecX q0_input = VecX(0),  // optional initial position
        const VecX q_d0_input = VecX(0) // optional initial velocity
    );

    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the NLP */
    bool get_nlp_info(
        Index&          n,
        Index&          m,
        Index&          nnz_jac_g,
        Index&          nnz_h_lag,
        IndexStyleEnum& index_style
    ) final override;

    /**@name Methods to block default compiler methods.
    *
    * The compiler automatically generates the following three methods.
    *  Since the default compiler implementation is generally not what
    *  you want (for all but the most simple classes), we usually
    *  put the declarations of these methods in the private section
    *  and never implement them. This prevents the compiler from
    *  implementing an incorrect "default" behavior without us
    *  knowing. (See Scott Meyers book, "Effective C++")
    */
    //@{
    DigitSingleStepOptimizer(
       const DigitSingleStepOptimizer&
    );

    DigitSingleStepOptimizer& operator=(
       const DigitSingleStepOptimizer&
    );

    std::shared_ptr<BezierCurves> bcPtr_;
    std::shared_ptr<Trajectories> trajPtr_; 

    std::shared_ptr<DigitConstrainedInverseDynamics> dcidPtr_;
    std::shared_ptr<ConstrainedInverseDynamics> cidPtr_;
    std::shared_ptr<InverseDynamics> idPtr_;
};

}; // namespace Digit
}; // namespace RAPTOR

#endif // DIGIT_SINGLE_STEP_OPTIMIZER_H
