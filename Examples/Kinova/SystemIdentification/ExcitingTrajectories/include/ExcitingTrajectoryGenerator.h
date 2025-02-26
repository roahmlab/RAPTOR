#ifndef CONDITION_NUMBER_OPTIMIZER_H
#define CONDITION_NUMBER_OPTIMIZER_H

#include "KinovaConstants.h"

#include "Optimizer.h"

#include "RegressorInverseDynamics.h"
#include "FixedFrequencyFourierCurves.h"

#include "TrajectoryTerminalConstraints.h"
#include "JointLimits.h"
#include "VelocityLimits.h"
#include "AccelerationLimits.h"
#include "TorqueLimits.h"
#include "KinovaCustomizedConstraints.h"

namespace RAPTOR {
namespace Kinova {

class ExcitingTrajectoryGenerator : public Optimizer {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    ExcitingTrajectoryGenerator() = default;

    /** Default destructor */
    ~ExcitingTrajectoryGenerator() = default;

    // [set_parameters]
    bool set_parameters(
        const VecX& x0_input,
        const Number T_input,
        const int N_input,
        const int degree_input,
        const double base_frequency_input,
        const VecX& q0_input,
        const VecX& q_d0_input,
        const Model& model_input, 
        const Eigen::VectorXi& independent_param_inds_input,
        const std::vector<Vec3>& boxCenters,
        const std::vector<Vec3>& boxOrientations,
        const std::vector<Vec3>& boxSizes,
        const VecX& joint_limits_buffer_input,
        const VecX& velocity_limits_buffer_input,
        const VecX& torque_limits_buffer_input,
        const bool include_gripper_or_not = false,
        const double collison_buffer_input = 0.0,
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
    ExcitingTrajectoryGenerator(
       const ExcitingTrajectoryGenerator&
    );

    ExcitingTrajectoryGenerator& operator=(
       const ExcitingTrajectoryGenerator&
    );

    Eigen::VectorXi independent_param_inds;

    std::shared_ptr<Trajectories> trajPtr_;

    std::shared_ptr<RegressorInverseDynamics> ridPtr_;

    MatX Y_independent;
    MatX pY_independent_pz;
};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // CONDITION_NUMBER_OPTIMIZER_H