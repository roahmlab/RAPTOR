#ifndef G1_SINGLE_STEP_OPTIMIZER_H
#define G1_SINGLE_STEP_OPTIMIZER_H

#include "Optimizer.h"

#include "BezierCurves.h"

#include "G1ConstrainedInverseDynamics.h"
#include "G1DynamicsConstraints.h"
#include "G1SingleStepPeriodicityConstraints.h"
#include "Utils.h"

#include "ConstrainedJointLimits.h"
#include "TorqueLimits.h"
#include "RectangleSurfaceContactConstraints.h"
#include "G1CustomizedConstraints.h"

#include "MinimizePower.h"
#include "MinimizeInitialVelocity.h"
#include "MinimizeInitialAcceleration.h"

namespace RAPTOR {
namespace G1 {

using namespace Ipopt;

class G1SingleStepOptimizer : public Optimizer {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    G1SingleStepOptimizer() = default;

    /** Default destructor */
    ~G1SingleStepOptimizer() = default;

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
        const Transform stance_foot_T_des = Transform(),
        const bool periodic = true,
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
    G1SingleStepOptimizer(
       const G1SingleStepOptimizer&
    );

    G1SingleStepOptimizer& operator=(
       const G1SingleStepOptimizer&
    );

    std::shared_ptr<BezierCurves> bcPtr_;
    std::shared_ptr<Trajectories> trajPtr_; 

    std::shared_ptr<DynamicsConstraints> dcPtr_;
    std::shared_ptr<InverseDynamics> idPtr_;
    std::shared_ptr<ConstrainedInverseDynamics> cidPtr_;
    std::shared_ptr<G1ConstrainedInverseDynamics> dcidPtr_;
};

}; // namespace G1
}; // namespace RAPTOR

#endif // G1_SINGLE_STEP_OPTIMIZER_H
