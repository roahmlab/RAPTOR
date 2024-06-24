#ifndef CONDITION_NUMBER_OPTIMIZER_H
#define CONDITION_NUMBER_OPTIMIZER_H

#include "KinovaConstants.h"

#include "Optimizer.h"

#include "FixedFrequencyFourierCurves.h"

#include "JointLimits.h"

namespace IDTO {
namespace Kinova {

class DataFilterOptimizer : public Optimizer {
public:
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    /** Default constructor */
    DataFilterOptimizer() = default;

    /** Default destructor */
    ~DataFilterOptimizer() = default;

    // [set_parameters]
    bool set_parameters(
        const VecX& x0_input,
        const VecX& tspan_input,
        const MatX& q_input,
        const MatX& q_d_input,
        const int degree_input,
        const int base_frequency_input
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

    /** Method to return the hessian of the objective */
    bool eval_hess_f(
        Index         n,
        const Number* x,
        bool          new_x,
        MatX&         hess_f
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
    DataFilterOptimizer(
       const DataFilterOptimizer&
    );

    DataFilterOptimizer& operator=(
       const DataFilterOptimizer&
    );

    const Number SQUARE_ROOT_THRESHOLD = 1e-8;

    MatX q_data;
    MatX q_d_data;

    Number position_weight = 5.0;
    Number velocity_weight = 1.0;

    std::shared_ptr<Trajectories> trajPtr_;
};

}; // namespace Kinova
}; // namespace IDTO

#endif // CONDITION_NUMBER_OPTIMIZER_H