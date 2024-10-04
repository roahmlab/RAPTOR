#ifndef FRICTIONPARAMETERSIDENTIFICATION_H
#define FRICTIONPARAMETERSIDENTIFICATION_H


#include "Optimizer.h"
#include "Trajectories.h" 

#include "RegressorInverseDynamics.h"
#include "FixedFrequencyFourierCurves.h"

#include "JointLimits.h"
#include "VelocityLimits.h"
#include "TorqueLimits.h"
#include "KinovaCustomizedConstraints.h"

namespace RAPTOR {
namespace Kinova {

class FrictionParametersIdentification : public Optimizer {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    FrictionParametersIdentification() = default;

    /** Default destructor */
    ~FrictionParametersIdentification() = default;

    // [set_parameters]
    bool set_parameters(
        VecX Xf,
        int nLinks,
        VecX& Fest,
        bool include_friction_offset,
        double N,
        std::shared_ptr<RegressorInverseDynamics>& RegressorID  
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

    bool get_bounds_info(
        Index n, 
        Number* x_l, 
        Number* x_u,
        Index m, 
        Number* g_l, 
        Number* g_u
    )final override;


    bool get_starting_point(
        Index       n, 
        bool        init_x, 
        Number*     x,
        bool        init_z, 
        Number*     z_L, 
        Number*     z_U,
        Index       m, 
        bool        init_lambda,
        Number*     lambda
    )final override;


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

    // bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) override;
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
    FrictionParametersIdentification(
       const FrictionParametersIdentification&
    );

    FrictionParametersIdentification& operator=(
       const FrictionParametersIdentification&
    );


    int nLinks_;
    VecX Fest_;
    bool include_friction_offset_;
    double N_;
    VecX Xf_;
    std::shared_ptr<RegressorInverseDynamics> ridPtr_;
   
};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // CONDITION_NUMBER_OPTIMIZER_H