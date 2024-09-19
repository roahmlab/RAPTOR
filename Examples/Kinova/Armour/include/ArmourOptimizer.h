#ifndef ARMOUR_OPTIMIZER_H
#define ARMOUR_OPTIMIZER_H

#include "KinovaConstants.h"

#include "ReachableSets.h"

namespace RAPTOR {
namespace Armour {
namespace Kinova {

class KinovaOptimizer : public Optimizer {
public:
    using Model = pinocchio::ModelTpl<float>;
    using VecX = Eigen::VectorXf;
    using Vec3 = Eigen::Vector3f;
    using MatX = Eigen::MatrixXf;

    /** Default constructor */
    KinovaOptimizer() = default;

    /** Default destructor */
    ~KinovaOptimizer() = default;

    // [set_parameters]
    bool set_parameters(
        const VecX& x0_input,
        const float T_input,
        const int N_input,
        const int degree_input,
        const Model& model_input, 
        const ArmourTrajectoryParameters& atp_input,
        const std::vector<Vec3>& boxCenters_input,
        const std::vector<Vec3>& boxOrientation_input,
        const std::vector<Vec3>& boxSize_input,
        const VecX& qdes_input,
        const int tplan_n_input,
        const VecX& joint_limits_buffer_input,
        const VecX& velocity_limits_buffer_input,
        const VecX& torque_limits_buffer_input,
        const float collision_buffer_input = 0
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
}; // namespace Armour
}; // namespace RAPTOR

#endif // ARMOUR_OPTIMIZER_H
