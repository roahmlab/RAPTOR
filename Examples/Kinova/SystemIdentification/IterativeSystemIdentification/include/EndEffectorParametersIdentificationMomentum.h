#ifndef ENDEFFECTOR_IDENTIFICATION_MOMENTUM_H
#define ENDEFFECTOR_IDENTIFICATION_MOMENTUM_H

#include "EndEffectorParametersIdentification.h"

namespace RAPTOR {

class EndEffectorParametersIdentificationMomentum : public EndEffectorParametersIdentification {
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
    EndEffectorParametersIdentificationMomentum() = default;

    /** Default destructor */
    ~EndEffectorParametersIdentificationMomentum() = default;

    // [set_parameters]
    bool set_parameters(
        const Model& model_input,
        const VecXd offset_input = VecXd::Zero(0)
    );

    // [add_trajectory_file]
    void add_trajectory_file(
        const std::string filename_input,
        const SensorNoiseInfo sensor_noise_input = SensorNoiseInfo(),
        const int H_input = 10,
        const TimeFormat time_format = TimeFormat::Second,
        const int downsample_rate = 1);

    // [initialize_regressors]
    void initialize_regressors(const std::shared_ptr<TrajectoryData>& trajPtr_,
                               const std::shared_ptr<TrajectoryData>& trajPtr2_,
                               const int H_input = 10);

    void reset() final override;

    /**@name Overloaded from TNLP */
    //@{

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
    EndEffectorParametersIdentificationMomentum(
       const EndEffectorParametersIdentificationMomentum&
    );

    EndEffectorParametersIdentificationMomentum& operator=(
       const EndEffectorParametersIdentificationMomentum&
    );

    // class members:
    std::vector<std::shared_ptr<TrajectoryData>> trajPtrs2_;

    std::shared_ptr<MomentumRegressor> mrPtr_;
    std::shared_ptr<RegressorInverseDynamics> ridPtr_;

    // std::shared_ptr<IntervalMomentumRegressor> mrIntPtr_ = nullptr;
    // std::shared_ptr<IntervalRegressorInverseDynamics> ridIntPtr_ = nullptr;

        // forward integration horizon
    int H = 100;

        // results
    Vec10d theta_uncertainty;
};

}; // namespace RAPTOR

#endif // ENDEFFECTOR_IDENTIFICATION_MOMENTUM_H