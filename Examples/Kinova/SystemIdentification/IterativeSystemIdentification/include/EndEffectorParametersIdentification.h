#ifndef ENDEFFECTOR_IDENTIFICATION_H
#define ENDEFFECTOR_IDENTIFICATION_H

#include "pinocchio/algorithm/regressor.hpp"

#include "Optimizer.h"
#include "MomentumRegressor.h"
#include "TrajectoryData.h"

namespace RAPTOR {

class EndEffectorParametersIdentification : public Optimizer {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Mat3 = Eigen::Matrix3d;
    using Mat4 = Eigen::Matrix4d;

    /** Default constructor */
    EndEffectorParametersIdentification() = default;

    /** Default destructor */
    ~EndEffectorParametersIdentification() = default;

    // [set_parameters]
    bool set_parameters(
        const Model& model_input,
        const std::string filename_input,
        const int H_input = 10,
        const VecX offset_input = VecX::Zero(0)
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
    EndEffectorParametersIdentification(
       const EndEffectorParametersIdentification&
    );

    EndEffectorParametersIdentification& operator=(
       const EndEffectorParametersIdentification&
    );

    // class members:
    std::shared_ptr<Model> modelPtr_; // robot model
    std::shared_ptr<Data> dataPtr_; // robot data

    VecX phi; // dynamic parameters of the robot model, the last 10 parameters are the end-effector parameters to be indentified
    VecX phi_original; // dynamic parameters read from the original robot model
   
    std::shared_ptr<TrajectoryData> trajPtr_;
    std::shared_ptr<TrajectoryData> trajPtr2_;

    std::shared_ptr<MomentumRegressor> mrPtr_;
    std::shared_ptr<RegressorInverseDynamics> ridPtr_;

        // forward integration horizon
    int H = 10;

        // offset in friction parameters
    bool include_offset = false;
    VecX offset;

        // regression data
    MatX A;
    VecX b;
};

}; // namespace RAPTOR

#endif // BASE_PARAMETERS_IDENTIFICATION_H