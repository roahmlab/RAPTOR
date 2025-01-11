#ifndef DUAL_ARMOUR_OPTIMIZER_H
#define DUAL_ARMOUR_OPTIMIZER_H

#include "ArmourOptimizer.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

class DualArmourOptimizer : public Optimizer {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;
    using Vec3 = Eigen::Vector3d;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    DualArmourOptimizer() = default;

    /** Default destructor */
    ~DualArmourOptimizer() = default;

    // [set_parameters]
    bool set_parameters(
        const VecX& q_des_input,
        Number t_plan_input,
        const std::shared_ptr<RobotInfo>& robotInfoPtr1_input,
        const std::shared_ptr<BezierCurveInterval>& trajPtr1_input,
        const std::shared_ptr<PZDynamics>& dynPtr1_input,
        const std::shared_ptr<RobotInfo>& robotInfoPtr2_input,
        const std::shared_ptr<BezierCurveInterval>& trajPtr2_input,
        const std::shared_ptr<PZDynamics>& dynPtr2_input,
        const std::vector<Vec3>& boxCenters_input,
        const std::vector<Vec3>& boxOrientation_input,
        const std::vector<Vec3>& boxSize_input
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
    DualArmourOptimizer(
       const DualArmourOptimizer&
    );

    DualArmourOptimizer& operator=(
       const DualArmourOptimizer&
    );

    std::shared_ptr<ArmourOptimizer> armourOptPtr1_ = nullptr;
    std::shared_ptr<ArmourOptimizer> armourOptPtr2_ = nullptr;

    std::shared_ptr<RobotInfo> robotInfoPtr1_ = nullptr;
    std::shared_ptr<RobotInfo> robotInfoPtr2_ = nullptr;

    std::vector<std::shared_ptr<TaperedCapsuleCollision<2 * NUM_FACTORS>>> tccPtrs;

    VecX g_lb_copy;
    VecX g_ub_copy;
};

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR

#endif // DUAL_ARMOUR_OPTIMIZER_H