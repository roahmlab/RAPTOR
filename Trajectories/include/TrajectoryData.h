#ifndef TRAJECTORY_DATA_H
#define TRAJECTORY_DATA_H

#include "Trajectories.h"

namespace RAPTOR {

// maximum noise on joint encoder sensor measurements
typedef struct SensorNoiseInfo_ {
    enum SensorNoiseType {
        Constant = 0, // constant noise
        Ratio // noise is a ratio of the measurement, for example 0.1 means 10% noise
    };

    SensorNoiseType position_error_type = SensorNoiseType::Constant;
    Eigen::VectorXd position_error;

    SensorNoiseType velocity_error_type = SensorNoiseType::Constant;
    Eigen::VectorXd velocity_error;

    SensorNoiseType acceleration_error_type = SensorNoiseType::Constant;
    Eigen::VectorXd acceleration_error;

    SensorNoiseInfo_() {
        position_error = Eigen::VectorXd::Zero(0);
        velocity_error = Eigen::VectorXd::Zero(0);
        acceleration_error = Eigen::VectorXd::Zero(0);
    };

    SensorNoiseInfo_(const int Nact) {
        position_error = Eigen::VectorXd::Zero(Nact);
        velocity_error = Eigen::VectorXd::Zero(Nact);
        acceleration_error = Eigen::VectorXd::Zero(Nact);
    };
} SensorNoiseInfo;

enum TimeFormat {
    Second = 0,
    Millisecond,
    Microsecond,
    Nanosecond
};

class TrajectoryData : public Trajectories {
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    TrajectoryData() = default;

    TrajectoryData(const std::string& filename_input,
                   const SensorNoiseInfo sensor_noise_input = SensorNoiseInfo(),
                   const TimeFormat time_format = TimeFormat::Second,
                   const int downsample_rate = 1,
                   const bool add_sensor_noise = false);

    TrajectoryData(const std::vector<std::string>& filenames_input,
                   const SensorNoiseInfo sensor_noise_input = SensorNoiseInfo(),
                   const TimeFormat time_format = TimeFormat::Second,
                   const int downsample_rate = 1,
                   const bool add_sensor_noise = false);

    TrajectoryData(double T_input,
                   int N_input, 
                   int Nact_input,
                   const bool random_or_not = true,
                   const SensorNoiseInfo sensor_noise_input = SensorNoiseInfo());

    // Destructor
    ~TrajectoryData() = default;

    // class methods:
        // has already loaded the data in the constructor
        // simply do nothing here
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) final override {};
                        
    // class members:
    SensorNoiseInfo sensor_noise;
};

}; // namespace RAPTOR

#endif // TRAJECTORY_DATA_H
