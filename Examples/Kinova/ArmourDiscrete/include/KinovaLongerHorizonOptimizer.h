#ifndef KINOVA_LONGER_HORIZON_OPTIMIZER_H
#define KINOVA_LONGER_HORIZON_OPTIMIZER_H

#include "KinovaConstants.h"

#include "Optimizer.h"

#include "PiecewiseBezierCurves.h"

#include "InverseDynamics.h"
#include "JointLimits.h"
#include "VelocityLimits.h"
#include "TorqueLimits.h"
#include "KinovaCustomizedConstraints.h"

#include "MinimizeTorque.h"
#include "MinimizePathLength.h"
#include "MinimizeJerk.h"

namespace RAPTOR {
namespace Kinova {

class KinovaLongerHorizonOptimizer : public Optimizer {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;
    using Vec3 = Eigen::Vector3d;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    KinovaLongerHorizonOptimizer() = default;

    /** Default destructor */
    ~KinovaLongerHorizonOptimizer() = default;

    // [set_parameters]
    bool set_parameters(
        const VecX& x0_input,
        const double T_input,
        const int N_input,
        const int degree_input,
        const Model& model_input,
        const VecX& q0_input,
        const VecX& qT_input,
        const std::vector<Vec3>& boxCenters_input,
        const std::vector<Vec3>& boxOrientation_input,
        const std::vector<Vec3>& boxSize_input,
        const VecX& joint_limits_buffer_input,
        const VecX& velocity_limits_buffer_input,
        const VecX& torque_limits_buffer_input,
        const bool include_gripper_or_not = false,
        const double collision_buffer_input = 0.0
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
    KinovaLongerHorizonOptimizer(
       const KinovaLongerHorizonOptimizer&
    );

    KinovaLongerHorizonOptimizer& operator=(
       const KinovaLongerHorizonOptimizer&
    );

    std::shared_ptr<Trajectories> trajPtr_;

    std::shared_ptr<InverseDynamics> idPtr_;

    VecX q_des;
    int tplan_n = 0;
};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // KINOVA_LONGER_HORIZON_OPTIMIZER_H
