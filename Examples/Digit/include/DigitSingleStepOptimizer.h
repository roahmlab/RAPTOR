#ifndef DIGITSINGLESTEPOPTIMIZER_H
#define DIGITSINGLESTEPOPTIMIZER_H

#include "Optimizer.h"

#include "FourierCurves.h"
#include "DigitConstrainedInverseDynamics.h"
#include "DigitDynamicsConstraints.h"
#include "Utils.h"

#include "ConstrainedJointLimits.h"
#include "TorqueLimits.h"
#include "SurfaceContactConstraints.h"
#include "KinematicsConstraints.h"

namespace IDTO {
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
        const int degree_input,
        const Model& model_input, 
        const Eigen::VectorXi& jtype_input,
        char stanceLeg, 
        const Transform& stance_foot_T_des_input
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

    /** Method to return the bounds for my problem */
    bool get_bounds_info(
        Index   n,
        Number* x_l,
        Number* x_u,
        Index   m,
        Number* g_l,
        Number* g_u
    ) override;

    /** Method to return the starting point for the algorithm */
    bool get_starting_point(
        Index   n,
        bool    init_x,
        Number* x,
        bool    init_z,
        Number* z_L,
        Number* z_U,
        Index   m,
        bool    init_lambda,
        Number* lambda
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

    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    void finalize_solution(
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
    ) override;
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
    DigitSingleStepOptimizer(
       const DigitSingleStepOptimizer&
    );

    DigitSingleStepOptimizer& operator=(
       const DigitSingleStepOptimizer&
    );

    std::shared_ptr<Trajectories> trajPtr_;
    std::shared_ptr<FourierCurves> fcPtr_;

    std::shared_ptr<DynamicsConstraints> dcPtr_;
    std::shared_ptr<InverseDynamics> idPtr_;
    std::shared_ptr<ConstrainedInverseDynamics> cidPtr_;
    std::shared_ptr<DigitConstrainedInverseDynamics> dcidPtr_;

    VecX x0;
};

}; // namespace Digit
}; // namespace IDTO

#endif // DIGITSINGLESTEPOPTIMIZER_H