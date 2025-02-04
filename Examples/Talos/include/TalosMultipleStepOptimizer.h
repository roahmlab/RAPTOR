#ifndef TALOS_MULTIPLE_STEP_OPTIMIZER_H
#define TALOS_MULTIPLE_STEP_OPTIMIZER_H

#include "TalosSingleStepOptimizer.h"
#include "TalosMultipleStepPeriodicityConstraints.h"

namespace RAPTOR {
namespace Talos {

using namespace Ipopt;

Eigen::VectorXd switchSolutionFromLeftToRight(const Eigen::VectorXd& z, 
                                              const int degree) {
    if (z.size() != (degree + 1) * NUM_INDEPENDENT_JOINTS + NUM_JOINTS + NUM_DEPENDENT_JOINTS) {
        throw std::invalid_argument("z has wrong size in switchSolutionFromLeftToRight! A single step solution is required.");
    }
    
    Eigen::VectorXd z_switched = z;

    // swap left leg and right leg
    z_switched.head((degree + 1) * NUM_INDEPENDENT_JOINTS / 2) = 
        z.segment((degree + 1) * NUM_INDEPENDENT_JOINTS / 2, (degree + 1) * NUM_INDEPENDENT_JOINTS / 2);
    z_switched.segment((degree + 1) * NUM_INDEPENDENT_JOINTS / 2, (degree + 1) * NUM_INDEPENDENT_JOINTS / 2) = 
        z.head((degree + 1) * NUM_INDEPENDENT_JOINTS / 2);

    return z_switched;
}

class TalosMultipleStepOptimizer : public Optimizer {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    TalosMultipleStepOptimizer() = default;

    /** Default destructor */
    ~TalosMultipleStepOptimizer() = default;

    // [set_parameters]
    bool set_parameters(
        const int NSteps_input,
        const VecX& x0_input,
        const double T_input,
        const int N_input,
        const TimeDiscretization time_discretization_input,
        const int degree_input,
        const Model& model_input, 
        const std::vector<GaitParameters>& gps_input
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
    TalosMultipleStepOptimizer(
       const TalosMultipleStepOptimizer&
    );

    TalosMultipleStepOptimizer& operator=(
       const TalosMultipleStepOptimizer&
    );

    std::vector<SmartPtr<TalosSingleStepOptimizer>> stepOptVec_;
    std::vector<std::shared_ptr<TalosMultipleStepPeriodicityConstraints>> periodConsVec_;
    std::vector<Index> n_local;
    std::vector<Index> n_position; // record the position of decision variables of each walking step in the full decision vector
    std::vector<Index> m_local;
    std::vector<Index> m_position; // record the position of constraints of each walking step in the full constraint vector
    bool ifFeasibleCurrIter = false;
};

}; // namespace Talos
}; // namespace RAPTOR

#endif // TALOS_MULTIPLE_STEP_OPTIMIZER_H
