#ifndef DIGIT_MULTIPLE_STEP_OPTIMIZER_H
#define DIGIT_MULTIPLE_STEP_OPTIMIZER_H

#include "DigitSingleStepOptimizer.h"
#include "DigitMultipleStepPeriodicityConstraints.h"

namespace RAPTOR {
namespace Digit {

using namespace Ipopt;

class DigitMultipleStepOptimizer : public Optimizer {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    DigitMultipleStepOptimizer() = default;

    /** Default destructor */
    ~DigitMultipleStepOptimizer() = default;

    // [set_parameters]
    bool set_parameters(
        const int NSteps_input,
        const VecX& x0_input,
        const double T_input,
        const int N_input,
        const TimeDiscretization time_discretization_input,
        const int degree_input,
        const Model& model_input, 
        const std::vector<GaitParameters>& gps_input,
        const bool periodic = false
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

    /** Method to return the bounds for my problem */
    bool get_bounds_info(
        Index   n,
        Number* x_l,
        Number* x_u,
        Index   m,
        Number* g_l,
        Number* g_u
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

    /** Method to return the constraint residuals */
    bool eval_g(
        Index         n,
        const Number* x,
        bool          new_x,
        Index         m,
        Number*       g
    ) final override;

    /** Method to return:
    *   1) The structure of the jacobian (if "values" is NULL)
    *   2) The values of the jacobian (if "values" is not NULL)
    */
    bool eval_jac_g(
        Index         n,
        const Number* x,
        bool          new_x,
        Index         m,
        Index         nele_jac,
        Index*        iRow,
        Index*        jCol,
        Number*       values
    ) final override;

    /** This method summarizes constraint violation for each type of constraints */
    void summarize_constraints(
        Index                      m,
        const Number*              g,
        const bool                 verbose = true
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
    DigitMultipleStepOptimizer(
       const DigitMultipleStepOptimizer&
    );

    DigitMultipleStepOptimizer& operator=(
       const DigitMultipleStepOptimizer&
    );

    // switch the solution from left stance to right stance
    // to be more specific, we have to perform the following changes:
    // 1. swap the left leg joints and right leg joints
    // 2. negate all joints
    // since Digit joint directions are mirrored between left leg and right leg
    VecX switchSolutionFromLeftToRight(std::shared_ptr<DigitConstrainedInverseDynamics>& currDcidPtr_,
                                       std::shared_ptr<DigitConstrainedInverseDynamics>& nextDcidPtr_,
                                       const VecX& z, 
                                       const int degree);

    std::vector<SmartPtr<DigitSingleStepOptimizer>> stepOptVec_;
    std::vector<std::shared_ptr<DigitMultipleStepPeriodicityConstraints>> periodConsVec_;
    std::vector<Index> n_local;
    std::vector<Index> n_position; // record the position of decision variables of each walking step in the full decision vector
    std::vector<Index> m_local;
    std::vector<Index> m_position; // record the position of constraints of each walking step in the full constraint vector
    bool ifFeasibleCurrIter = false;
};

}; // namespace Digit
}; // namespace RAPTOR

#endif // DIGIT_MULTIPLE_STEP_OPTIMIZER_H
