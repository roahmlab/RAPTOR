#ifndef KINOVA_IK_SOLVER_H
#define KINOVA_IK_SOLVER_H

#include "KinovaConstants.h"

#include "Optimizer.h"

#include "Plain.h"

#include "JointLimits.h"
#include "KinematicsConstraints.h"
#include "KinovaCustomizedConstraints.h"

namespace RAPTOR {
namespace Kinova {

class KinovaIKSolver : public Optimizer {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;
    using Vec3 = Eigen::Vector3d;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    KinovaIKSolver() = default;

    /** Default destructor */
    ~KinovaIKSolver() = default;

    // [set_parameters]
    bool set_parameters(
        const VecX& x0_input,
        const Model& model_input,
        const Transform& desiredTransform_input,
        const std::vector<Vec3>& boxCenters_input,
        const std::vector<Vec3>& boxOrientation_input,
        const std::vector<Vec3>& boxSize_input,
        const Transform endT_input = Transform(),
        const bool include_gripper_or_not = true,
        const double collision_buffer_input = 0,
        Eigen::VectorXi jtype_input = Eigen::VectorXi(0)
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

    /** Method to return the hessian of the objective */
    bool eval_hess_f(
        Index         n,
        const Number* x,
        bool          new_x,
        MatX&         hess_f
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
    KinovaIKSolver(
       const KinovaIKSolver&
    );

    KinovaIKSolver& operator=(
       const KinovaIKSolver&
    );

    std::shared_ptr<Trajectories> trajPtr_;
    VecX start;

};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // KINOVA_IK_SOLVER_H
