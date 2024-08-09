#ifndef TALOSSINGLESTEPOPTIMIZER_H
#define TALOSSINGLESTEPOPTIMIZER_H

#include "Optimizer.h"

#include "BezierCurves.h"

#include "TalosConstrainedInverseDynamics.h"
#include "TalosDynamicsConstraints.h"
#include "Utils.h"

#include "ConstrainedJointLimits.h"
#include "TorqueLimits.h"
#include "SurfaceContactConstraints.h"
#include "TalosCustomizedConstraints.h"

namespace RAPTOR {
namespace Talos {

using namespace Ipopt;

class TalosSingleStepOptimizer : public Optimizer {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    TalosSingleStepOptimizer() = default;

    /** Default destructor */
    ~TalosSingleStepOptimizer() = default;

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

    /** Method to return the objective value */
    bool eval_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number&       obj_value
    ) final override;

    /** Method to return the gradient of the objective */
    bool eval_grad_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number*       grad_f
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
    TalosSingleStepOptimizer(
       const TalosSingleStepOptimizer&
    );

    TalosSingleStepOptimizer& operator=(
       const TalosSingleStepOptimizer&
    );

    std::shared_ptr<BezierCurves> bcPtr_;
    std::shared_ptr<Trajectories> trajPtr_; 

    std::shared_ptr<DynamicsConstraints> dcPtr_;
    std::shared_ptr<InverseDynamics> idPtr_;
    std::shared_ptr<ConstrainedInverseDynamics> cidPtr_;
    std::shared_ptr<TalosConstrainedInverseDynamics> dcidPtr_;
};

}; // namespace Talos
}; // namespace RAPTOR

#endif // TALOSSINGLESTEPOPTIMIZER_H
