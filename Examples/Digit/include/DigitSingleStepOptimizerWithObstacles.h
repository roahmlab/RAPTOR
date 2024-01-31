#ifndef DIGITSINGLESTEPOPTIMIZER_WITHOBSTACLES_H
#define DIGITSINGLESTEPOPTIMIZER_WITHOBSTACLES_H

#include "Optimizer.h"

#include "FourierCurves.h"
#include "FixedFrequencyFourierCurves.h"
#include "BezierCurves.h"
#include "DigitConstrainedInverseDynamics.h"
#include "DigitDynamicsConstraints.h"
#include "Utils.h"

#include "ConstrainedJointLimits.h"
#include "TorqueLimits.h"
#include "SurfaceContactConstraints.h"
#include "DigitCustomizedConstraintsWithObstacles.h"

namespace IDTO {
namespace Digit {

using namespace Ipopt;

class DigitSingleStepOptimizerWithObstacles : public Optimizer {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    DigitSingleStepOptimizerWithObstacles() = default;

    /** Default destructor */
    ~DigitSingleStepOptimizerWithObstacles() = default;

    // [set_parameters]
    bool set_parameters(
        const VecX& x0_input,
        const double T_input,
        const int N_input,
        const TimeDiscretization time_discretization_input,
        const int degree_input,
        const Model& model_input, 
        const Eigen::VectorXi& jtype_input,
        const Eigen::Array<Vec3, 1, Eigen::Dynamic>& zonotopeCenters_input,
        const Eigen::Array<MatX, 1, Eigen::Dynamic>& zonotopeGenerators_input,
        const VecX& q_act0_input,
        const VecX& q_act_d0_input 
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
    ) override;

    /** Method to return the objective value */
    bool eval_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number&       obj_value
    ) override;

    /** Method to return the gradient of the objective */
    bool eval_grad_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number*       grad_f
    ) override;

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
    DigitSingleStepOptimizerWithObstacles(
       const DigitSingleStepOptimizerWithObstacles&
    );

    DigitSingleStepOptimizerWithObstacles& operator=(
       const DigitSingleStepOptimizerWithObstacles&
    );

    std::shared_ptr<Trajectories> trajPtr_; 

    std::shared_ptr<ConstrainedInverseDynamics> cidPtr_;
};

}; // namespace Digit
}; // namespace IDTO

#endif // DIGITSINGLESTEPOPTIMIZER_WITHOBSTACLES_H
