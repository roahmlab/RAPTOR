#ifndef BASEPARAMETERSIDENTIFICATION_H
#define BASEPARAMETERSIDENTIFICATION_H


#include "Optimizer.h"

#include "RegressorInverseDynamics.h"
#include "FixedFrequencyFourierCurves.h"

#include "JointLimits.h"
#include "VelocityLimits.h"
#include "TorqueLimits.h"
#include "KinovaCustomizedConstraints.h"
#include "QRDecompositionSolver.h"

namespace RAPTOR {
namespace Kinova {

class BaseParametersIdentification : public Optimizer {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Mat3 = Eigen::Matrix3d;


    /** Default constructor */
    BaseParametersIdentification() = default;

    /** Default destructor */
    ~BaseParametersIdentification() = default;

    // [set_parameters]
    bool set_parameters(
        MatX &Wh,
        VecX &Th,
        VecX &X,
        bool include_friction_offset,
        Model &model_input,
        std::shared_ptr<QRDecompositionSolver> regroupPtr,
        VecX &lb,
        VecX &ub,
        int b_full,
        int fm_dim,
        int Alg_case
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

    
    bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                                   Index m, Number* g_l, Number* g_u)final override;

    bool get_starting_point(Index n, bool init_x, Number* x,
                                              bool init_z, Number* z_L, Number* z_U,
                                              Index m, bool init_lambda,
                                              Number* lambda)final override;

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

    bool eval_g(Index n, const Number *x, bool new_x,
                                                Index m, Number *g)final override;

    bool eval_jac_g(Index n, const Number *x, bool new_x,
                    Index m, Index nele_jac, Index *iRow, Index *jCol,
                    Number *values)final override;

    // bool eval_h(
    //     Index n, const Number* x, bool new_x, Number obj_factor,
    //     Index m, const Number* lambda, bool new_lambda,
    //     Index nele_hess, Index* iRow, Index* jCol, Number* values)final override;

    void compute_LMI_matrix(const VecX &pi_inertia,  Index j, MatX &LMI);
    void compute_LMI_gradient(const VecX &pi_full, Index j, MatX &dLMIdpi_full);




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
    BaseParametersIdentification(
       const BaseParametersIdentification&
    );

    BaseParametersIdentification& operator=(
       const BaseParametersIdentification&
    );


    MatX Wh_;
    VecX Th_;
    VecX X_;
    int fm_dim_;
    bool include_friction_offset_;
    Model model_input_;
    VecX lb_;
    VecX ub_;
    int b_full_;
    int Alg_case_;
    int nLinks_;

    // regroup 
    MatX Ginv_;
    VecX pi_d_;
    int b_dim_;
    int d_dim_;
    int p_ip_;

    std::shared_ptr<QRDecompositionSolver> regroupPtr_;
};

} // namespace Kinova
} // namespace RAPTOR

#endif // BASEPARAMETERSIDENTIFICATION_H