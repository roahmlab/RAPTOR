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
    using VecX = Eigen::VectorXd;
    using Vec10 = Eigen::Vector<double, 10>;
    using MatX = Eigen::MatrixXd;
    using Mat4 = Eigen::Matrix4d;
    using Mat10 = Eigen::Matrix<double, 10, 10>;

    /** Default constructor */
    EndEffectorParametersIdentification() = default;

    /** Default destructor */
    ~EndEffectorParametersIdentification() = default;

    // [set_parameters]
    bool set_parameters(
        const Model& model_input,
        const VecX offset_input = VecX::Zero(0)
    );

    // [add_trajectory_file]
    void add_trajectory_file(
        const std::shared_ptr<TrajectoryData>& trajectory_input,
        const MatX& acceleration_input);

    void add_trajectory_file(
        const std::string hardware_trajectory_filename_input,
        const std::string acceleration_filename_input,
        const TimeFormat time_format = TimeFormat::Second,
        const int downsample_rate = 1);

    // [initialize_regressors]
    void initialize_regressors(const std::shared_ptr<TrajectoryData>& trajPtr_,
                               const MatX& acceleration);

    void reset() override;

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
    Vec10 z_to_theta(const VecX& z);

    Vec10 d_z_to_theta(
        const VecX& z,
        Mat10& dtheta);

    Vec10 dd_z_to_theta(
        const VecX& z,
        Mat10& dtheta,
        Eigen::Array<Mat10, 1, 10>& ddtheta);

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
   
    std::vector<std::shared_ptr<TrajectoryData>> trajPtrs_;
    std::vector<MatX> accelerations_;
    std::vector<std::string> trajectoryFilenames_;

        // number of segments / trajectory data files
    std::vector<int> num_segments;
    int total_num_segments = 0;

        // offset in friction parameters
    bool include_offset = false;
    VecX offset;

        // regression data
    std::vector<MatX> Aseg; // regression matrix for each trajectory
    std::vector<VecX> bseg; // regression vector for each trajectory

        // results
    Vec10 theta_solution;
};

}; // namespace RAPTOR

#endif // ENDEFFECTOR_IDENTIFICATION_H