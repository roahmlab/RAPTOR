#ifndef KINOVA_OPTIMIZER_H
#define KINOVA_OPTIMIZER_H

#include "KinovaCustomizedConstraints.h"
#include "KinovaConstants.h"
#include "InverseDynamics.h"
#include "Optimizer.h"
#include "ArmourBezierCurves.h"
#include "JointLimits.h"
#include "VelocityLimits.h"
#include "TorqueLimits.h"

namespace IDTO {
namespace Kinova {

class KinovaOptimizer : public Optimizer {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;
    using Vec3 = Eigen::Vector3d;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    KinovaOptimizer() = default;

    /** Default destructor */
    ~KinovaOptimizer() = default;

    // [set_parameters]
    bool set_parameters(
        const VecX& x0_input,
        const double T_input,
        const int N_input,
        const int degree_input,
        const Model& model_input, 
        const Eigen::VectorXi& jtype_input,
        const ArmourTrajectoryParameters& atp_input,
        const std::vector<Vec3>& boxCenters_input,
        const std::vector<Vec3>& boxOrientation_input,
        const std::vector<Vec3>& boxSize_input,
        const VecX& qdes_input,
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
    KinovaOptimizer(
       const KinovaOptimizer&
    );

    KinovaOptimizer& operator=(
       const KinovaOptimizer&
    );

    std::shared_ptr<Trajectories> trajPtr_;

    std::shared_ptr<InverseDynamics> idPtr_;

    VecX qdes;
    int tplan_n = 0;
};

}; // namespace Kinova
}; // namespace IDTO

#endif // KINOVA_OPTIMIZER_H
