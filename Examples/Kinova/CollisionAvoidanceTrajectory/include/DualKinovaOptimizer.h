#ifndef DUAL_KINOVA_OPTIMIZER_H
#define DUAL_KINOVA_OPTIMIZER_H

#include "KinovaLongerHorizonOptimizer.h"

#include "TaperedCapsuleCollision.h"

namespace RAPTOR {
namespace Kinova {

class DualKinovaOptimizer : public Optimizer {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;
    using Vec3 = Eigen::Vector3d;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    DualKinovaOptimizer() = default;

    /** Default destructor */
    ~DualKinovaOptimizer() = default;

    // [set_parameters]
    bool set_parameters(
        const VecX& x0_input,
        const double T_input,
        const int N_input,
        const Model& model1_input, 
        const Model& model2_input, 
        const VecX& q0_input,
        const VecX& qT_input,
        const std::vector<Vec3>& boxCenters_input,
        const std::vector<Vec3>& boxOrientation_input,
        const std::vector<Vec3>& boxSize_input,
        const VecX& joint_limits_buffer_input,
        const VecX& velocity_limits_buffer_input,
        const VecX& torque_limits_buffer_input,
        const bool include_gripper_or_not = true,
        const double collision_buffer_input = 0.0,
        const double arm_arm_collision_buffer_input = 0.0
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

    void summarize_constraints(
        Index                      m,
        const Number*              g,
        const bool                 verbose
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
    DualKinovaOptimizer(
       const DualKinovaOptimizer&
    );

    DualKinovaOptimizer& operator=(
       const DualKinovaOptimizer&
    );

    std::shared_ptr<KinovaLongerHorizonOptimizer> kinovaOptPtr1_;
    std::shared_ptr<KinovaLongerHorizonOptimizer> kinovaOptPtr2_;

    const static int degree_input = 1;
    const static int numVars_dummy = 7 * degree_input * 3 * 2;

    TaperedCapsuleCollision<numVars_dummy> tcc;
    double arm_arm_collision_buffer = 0.0;

    VecX g_lb_copy;
    VecX g_ub_copy;
};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // DUAL_KINOVA_OPTIMIZER_H
