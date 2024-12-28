#ifndef ENDEFFECTOR_IDENTIFICATION_H
#define ENDEFFECTOR_IDENTIFICATION_H

#include "pinocchio/algorithm/regressor.hpp"

#include "Optimizer.h"
#include "MomentumRegressor.h"
#include "IntervalMomentumRegressor.h"
#include "TrajectoryData.h"

namespace RAPTOR {

class EndEffectorParametersIdentification : public Optimizer {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecXd = Eigen::VectorXd;
    using Vec10d = Eigen::Vector<double, 10>;
    using MatXd = Eigen::MatrixXd;
    using Mat4d = Eigen::Matrix4d;
    using Mat10d = Eigen::Matrix<double, 10, 10>;
    using VecXInt = Eigen::Vector<Interval, Eigen::Dynamic>;
    using MatXInt = Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic>;

    /** Default constructor */
    EndEffectorParametersIdentification() = default;

    /** Default destructor */
    ~EndEffectorParametersIdentification() = default;

    // [set_parameters]
    bool set_parameters(
        const Model& model_input,
        const std::string filename_input,
        const SensorNoiseInfo sensor_noise_input = SensorNoiseInfo(),
        const int H_input = 10,
        const TimeFormat time_format = TimeFormat::Second,
        const int downsample_rate = 1,
        const VecXd offset_input = VecXd::Zero(0)
    );

    // [set_parameters]
    bool set_parameters(
        const Model& model_input,
        const std::vector<std::string>& filenames_input,
        const SensorNoiseInfo sensor_noise_input = SensorNoiseInfo(),
        const int H_input = 10,
        const TimeFormat time_format = TimeFormat::Second,
        const int downsample_rate = 1,
        const VecXd offset_input = VecXd::Zero(0)
    );

    // [initialize_regressors]
    void initialize_regressors(const std::shared_ptr<TrajectoryData>& trajPtr_,
                               const std::shared_ptr<TrajectoryData>& trajPtr2_,
                               const int H_input = 10);

    void reset() final override;

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

    /** convert the decision variable to the dynamic parameters of the end effector */
    Vec10d z_to_theta(const VecXd& z);

    Vec10d d_z_to_theta(
        const VecXd& z,
        Mat10d& dtheta);

    Vec10d dd_z_to_theta(
        const VecXd& z,
        Mat10d& dtheta,
        Eigen::Array<Mat10d, 1, 10>& ddtheta);

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

    void finalize_solution(
        SolverReturn               status,
        Index                      n,
        const Number*              x,
        const Number*              z_L,
        const Number*              z_U,
        Index                      m,
        const Number*              g,
        const Number*              lambda,
        Number                     obj_value,
        const IpoptData*           ip_data,
        IpoptCalculatedQuantities* ip_cq
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

    VecXd phi; // dynamic parameters of the robot model, the last 10 parameters are the end-effector parameters to be indentified
    VecXd phi_original; // dynamic parameters read from the original robot model
   
    std::vector<std::shared_ptr<TrajectoryData>> trajPtrs_;
    std::vector<std::shared_ptr<TrajectoryData>> trajPtrs2_;

    std::shared_ptr<MomentumRegressor> mrPtr_;
    std::shared_ptr<RegressorInverseDynamics> ridPtr_;

    // std::shared_ptr<IntervalMomentumRegressor> mrIntPtr_ = nullptr;
    // std::shared_ptr<IntervalRegressorInverseDynamics> ridIntPtr_ = nullptr;

    VecXd weights; // weights for the nonlinear least square problem

        // forward integration horizon
    int H = 10;
    std::vector<int> num_segments;

        // offset in friction parameters
    bool include_offset = false;
    VecXd offset;

        // regression data
    std::vector<MatXd> Aseg; // regression matrix for each trajectory
    std::vector<VecXd> bseg; // regression vector for each trajectory

    MatXd A; // regression matrix for all trajectories
    VecXd b; // regression vector for all trajectories

    MatXd Aweighted;
    VecXd bweighted;

    Index nonzero_weights = 0;

        // results
    Vec10d theta_solution;
    Vec10d theta_uncertainty;
};

}; // namespace RAPTOR

#endif // BASE_PARAMETERS_IDENTIFICATION_H