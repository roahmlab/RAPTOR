#ifndef KINOVA_WAITR_OPTIMIZER_H
#define KINOVA_WAITR_OPTIMIZER_H

#include "pinocchio/algorithm/model.hpp"

#include "KinovaConstants.h"

#include "ArmourBezierCurves.h"

#include "Optimizer.h"

#include "JointLimits.h"
#include "VelocityLimits.h"
#include "TorqueLimits.h"
#include "CustomizedInverseDynamics.h"
#include "KinovaCustomizedConstraints.h"
#include "CircleSurfaceContactConstraints.h"

namespace RAPTOR {
namespace Kinova {

class KinovaWaitrOptimizer : public Optimizer {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    KinovaWaitrOptimizer() = default;

    /** Default destructor */
    ~KinovaWaitrOptimizer() = default;

    // [set_parameters]
    bool set_parameters(
        const VecX& x0_input,
        const double T_input,
        const int N_input,
        const int degree_input,
        const Model& model_input, 
        const ArmourTrajectoryParameters& atp_input,
        const circleContactSurfaceParams& csp_input,
        const std::vector<Vec3>& boxCenters_input,
        const std::vector<Vec3>& boxOrientation_input,
        const std::vector<Vec3>& boxSize_input,
        const VecX& q_des_input,
        const int tplan_n_input,
        const VecX& joint_limits_buffer_input,
        const VecX& velocity_limits_buffer_input,
        const VecX& torque_limits_buffer_input
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
    KinovaWaitrOptimizer(
       const KinovaWaitrOptimizer&
    );

    KinovaWaitrOptimizer& operator=(
       const KinovaWaitrOptimizer&
    );

    std::shared_ptr<Trajectories> trajPtr_;

    std::shared_ptr<CustomizedInverseDynamics> idPtr_;

    VecX q_des;
    int tplan_n = 0;
};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // KINOVA_WAITR_OPTIMIZER_H
