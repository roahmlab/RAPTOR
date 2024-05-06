#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"

#include <memory>
#include <map>
#include <string>
#include <chrono>
#include <omp.h>

#include "Constraints.h"

namespace IDTO {

using namespace Ipopt;

class Optimizer : public Ipopt::TNLP {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd; 

    /** Default constructor */
    Optimizer() = default;

    /** Default destructor */
    virtual ~Optimizer() = default;

    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the NLP */
    virtual bool get_nlp_info(
        Index&          n,
        Index&          m,
        Index&          nnz_jac_g,
        Index&          nnz_h_lag,
        IndexStyleEnum& index_style
    ) = 0;

    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(
        Index   n,
        Number* x_l,
        Number* x_u,
        Index   m,
        Number* g_l,
        Number* g_u
    );

    /** Method to return the starting point for the algorithm */
    virtual bool get_starting_point(
        Index   n,
        bool    init_x,
        Number* x,
        bool    init_z,
        Number* z_L,
        Number* z_U,
        Index   m,
        bool    init_lambda,
        Number* lambda
    );

    /** Method to warm start the entire iterate */
    virtual bool GetWarmStartIterate(
        IteratesVector& /*warm_start_iterate*/
    )
    {
        return false;
    }

    /** Method to return the objective value */
    virtual bool eval_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number&       obj_value
    ) = 0;

    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(
        Index         n,
        const Number* x,
        bool          new_x,
        Number*       grad_f
    ) = 0;

    /** Method to return the constraint residuals */
    virtual bool eval_g(
        Index         n,
        const Number* x,
        bool          new_x,
        Index         m,
        Number*       g
    );

    /** Method to return:
    *   1) The structure of the jacobian (if "values" is NULL)
    *   2) The values of the jacobian (if "values" is not NULL)
    */
    virtual bool eval_jac_g(
        Index         n,
        const Number* x,
        bool          new_x,
        Index         m,
        Index         nele_jac,
        Index*        iRow,
        Index*        jCol,
        Number*       values
    );

    /** Method to return:
    *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
    *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
    */
    virtual bool eval_h(
        Index         n,
        const Number* x,
        bool          new_x,
        Number        obj_factor,
        Index         m,
        const Number* lambda,
        bool          new_lambda,
        Index         nele_hess,
        Index*        iRow,
        Index*        jCol,
        Number*       values
    );

    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(
        SolverReturn               status,
        Index                      n,
        const Number*              x,
        const Number*              z_L,
        const Number*              z_U,
        Index                      m,
        const Number*              g,
        const Number*              lambda,
        Number                     obj_value,
        const IpoptData*           ip_data,
        IpoptCalculatedQuantities* ip_cq
    );
    //@}

    /** This method summarizes constraint violation for each type of constraints */
    virtual void summarize_constraints(
        Index                      m,
        const Number*              g
    );
    //@}

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
    Optimizer(
       const Optimizer&
    );

    Optimizer& operator=(
       const Optimizer&
    );

    // class members:
    int numVars = 0; // number of variables
    int numCons = 0; // number of constraints

    VecX x0; // stores the initial guess here

    std::vector<std::unique_ptr<Constraints>> constraintsPtrVec_;
    std::vector<std::string> constraintsNameVec_;

    VecX solution; // stores the final solution here

    Number final_constr_violation = 0;
    bool ifFeasible = true;
};

}; // namespace IDTO

#endif // OPTIMIZER_H