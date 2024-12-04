#ifndef ENDEFFECTOR_MOMENT_H
#define ENDEFFECTOR_MOMENT_H

#include "Optimizer.h"
// #include "QRDecompositionSolver.h"
#include "LMIConstraints.h"

#include "pinocchio/algorithm/regressor.hpp"

namespace RAPTOR {

class EndEffectorMoment : public Optimizer {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Mat3 = Eigen::Matrix3d;
    using Mat4 = Eigen::Matrix4d;
    

    /** Default constructor */
    EndEffectorMoment() = default;

    /** Default destructor */
    ~EndEffectorMoment() = default;

    // [set_parameters]
    bool set_parameters(
        const std::shared_ptr<MatX>& A_full_input,
        const std::shared_ptr<VecX>& b_input,
        const std::shared_ptr<VecX>& inertia_parametersPtr_input,
        const bool include_offset_input
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

    // /** Method to return the bounds for my problem */
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
    EndEffectorMoment(
       const EndEffectorMoment&
    );

    EndEffectorMoment& operator=(
       const EndEffectorMoment&
    );


    // shared pointers to data
    // std::shared_ptr<MatX> A_;
    // std::shared_ptr<VecX> b_;
    std::shared_ptr<VecX> inertia_parametersPtr_;

    MatX A;
    VecX b;
   
    MatX tau_inertials; // computed from the trajectory without friction
    int Nact = 0; // number of motors
    int N = 0; // number of samples

    bool include_offset = false;
};

}; // namespace RAPTOR

#endif // BASE_PARAMETERS_IDENTIFICATION_H