#ifndef ARMOUR_OPTIMIZER_H
#define ARMOUR_OPTIMIZER_H

#include "PZDynamics.h"
#include "BoxCollisionAvoidance.h"
#include "Optimizer.h"

#include "ArmourBezierCurves.h"
#include "EndEffectorRegressorConditionNumber.h"
#include "MomentumRegressor.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

class SafePayloadExcitingTrajectoryGenerator : public Optimizer {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;
    using Vec3 = Eigen::Vector3d;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    SafePayloadExcitingTrajectoryGenerator() = default;

    /** Default destructor */
    ~SafePayloadExcitingTrajectoryGenerator() = default;

    // [set_parameters]
    bool set_parameters(
        const std::shared_ptr<RobotInfo>& robotInfoPtr_input,
        const std::shared_ptr<BezierCurveInterval>& trajPtr_input,
        const std::shared_ptr<PZDynamics>& dynPtr_input,
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
    SafePayloadExcitingTrajectoryGenerator(
       const SafePayloadExcitingTrajectoryGenerator&
    );

    SafePayloadExcitingTrajectoryGenerator& operator=(
       const SafePayloadExcitingTrajectoryGenerator&
    );

    std::shared_ptr<RobotInfo> robotInfoPtr_;
    std::shared_ptr<BezierCurveInterval> trajIntervalPtr_;
    std::shared_ptr<PZDynamics> dynPtr_;
    
    std::vector<std::shared_ptr<BoxCollisionAvoidance>> bcaPtrs;

    size_t num_time_steps = 0;
    size_t num_spheres = 0;
    size_t num_fixed_joints = 0;

    std::shared_ptr<Trajectories> trajPtr_;
    std::shared_ptr<RegressorInverseDynamics> ridPtr_;
    std::shared_ptr<Costs> costPtr_;

    VecX g_lb_copy;
    VecX g_ub_copy;
};

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR

#endif // ARMOUR_OPTIMIZER_H
