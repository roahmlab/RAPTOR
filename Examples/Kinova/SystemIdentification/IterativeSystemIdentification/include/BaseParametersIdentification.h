#ifndef BASE_PARAMETERS_IDENTIFICATION_H
#define BASE_PARAMETERS_IDENTIFICATION_H

#include "Optimizer.h"
#include "QRDecompositionSolver.h"
#include "RegroupedLMIConstraints.h"

#include "pinocchio/algorithm/regressor.hpp"

namespace RAPTOR {

class BaseParametersIdentification : public Optimizer {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
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
        const Model& model_input,
        const std::shared_ptr<MatX>& posPtr_input,
        const std::shared_ptr<MatX>& velPtr_input,
        const std::shared_ptr<MatX>& accPtr_input,
        const std::shared_ptr<MatX>& torquePtr_input,
        std::shared_ptr<QRDecompositionSolver> regroupPtr_input,
        const bool include_offset_input = false
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
    BaseParametersIdentification(
       const BaseParametersIdentification&
    );

    BaseParametersIdentification& operator=(
       const BaseParametersIdentification&
    );

    const double default_maximum_uncertainty = 1; // default maximum uncertainty

    std::shared_ptr<Model> modelPtr_; // robot model
    std::shared_ptr<Data> dataPtr_; // robot data

    MatX FullObservationMatrix; // full observation matrix
    MatX RegroupedObservationMatrix; // regrouped observation matrix

    // shared pointers to data
    std::shared_ptr<MatX> posPtr_;
    std::shared_ptr<MatX> velPtr_;
    std::shared_ptr<MatX> accPtr_;
    std::shared_ptr<MatX> torquePtr_;

    MatX tau_inertials; // computed from the trajectory without friction

    int Nact = 0; // number of motors
    int N = 0; // number of samples

    std::shared_ptr<QRDecompositionSolver> regroupPtr_;

    bool include_offset = false;
};

}; // namespace RAPTOR

#endif // BASE_PARAMETERS_IDENTIFICATION_H