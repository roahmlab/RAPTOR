#include "TrajectoryData.h"

namespace RAPTOR {

TrajectoryData::TrajectoryData(const std::string& filename_input,
                               const SensorNoiseInfo sensor_noise_input,
                               const int downsample_rate) :
    sensor_noise(sensor_noise_input) {
    // parse file
    MatX traj_data = Utils::initializeEigenMatrixFromFile(filename_input);

    // check downsample format
    if (downsample_rate <= 0) {
        throw std::invalid_argument("Invalid downsample rate");
    }
    else if (downsample_rate > 1) {
        int num_samples = traj_data.rows() / downsample_rate;
        std::cout << "Performing downsample from " << traj_data.rows() << " to " << num_samples << std::endl;
        MatX new_traj_data(num_samples, traj_data.cols());
        Utils::uniformlySampleMatrixInRows(traj_data, new_traj_data, num_samples);
        traj_data = new_traj_data; // this operation can not be reversed!
    }

    N = traj_data.rows();
    Nact = (traj_data.cols() - 1) / 3;

    if (Nact <= 0) {
        throw std::invalid_argument("0 actuated joints");
    }
   
    if (2 * Nact + 1 != traj_data.cols() &&
        3 * Nact + 1 != traj_data.cols()) {
        std::cerr << Nact << ' ' << traj_data.cols() << std::endl;
        throw std::invalid_argument("Invalid trajectory file format");
    }

    tspan = traj_data.col(0);
    T = tspan(N - 1);

    // output information
    std::cout << "TrajectoryData: " << N << " data points loaded from " << filename_input << std::endl;
    std::cout << "TrajectoryData: T = " << T << std::endl;
    std::cout << "TrajectoryData: Nact = " << Nact << std::endl;

    // no variable length parameters, we don't compute gradient here
    // trajectory data is usually large so save some memory here
    varLength = 0;
    initialize_memory();

    for (int i = 0; i < N; i++) {
        q(i) = traj_data.row(i).segment(1, Nact);
        q_d(i) = traj_data.row(i).segment(1 + Nact, Nact);
    }

    // this could be acceleration estimation or applied torque from the sensor
    if (3 * Nact + 1 == traj_data.cols()) {
        std::cout << "TrajectoryData: q_dd or torque is loaded from the file" << std::endl;
        for (int i = 0; i < N; i++) {
            q_dd(i) = traj_data.row(i).segment(1 + 2 * Nact, Nact);
        }
    }

    // check sensor noise format
    if (sensor_noise.position_error.size() == 0 &&
        sensor_noise.velocity_error.size() == 0 &&
        sensor_noise.acceleration_error.size() == 0) {
        sensor_noise.position_error = VecX::Zero(Nact);
        sensor_noise.velocity_error = VecX::Zero(Nact);
        sensor_noise.acceleration_error = VecX::Zero(Nact);    
    }
    else if (sensor_noise.position_error.size() != Nact ||
        sensor_noise.velocity_error.size() != Nact ||
        sensor_noise.acceleration_error.size() != Nact) {
        throw std::invalid_argument("Invalid sensor noise format");
    }
}

TrajectoryData::TrajectoryData(double T_input,
                               int N_input, 
                               int Nact_input,
                               const SensorNoiseInfo sensor_noise_input) :
    Trajectories(0, T_input, N_input, Nact_input, TimeDiscretization::Uniform),
    sensor_noise(sensor_noise_input) {
    // randomly generate trajectory data
    for (int i = 0; i < N; i++) {
        q(i).setRandom();
        q_d(i).setRandom();
        q_dd(i).setRandom();
    }

    // check sensor noise format
    if (sensor_noise.position_error.size() == 0 &&
        sensor_noise.velocity_error.size() == 0 &&
        sensor_noise.acceleration_error.size() == 0) {
        sensor_noise.position_error = VecX::Zero(Nact);
        sensor_noise.velocity_error = VecX::Zero(Nact);
        sensor_noise.acceleration_error = VecX::Zero(Nact);    
    }
    else if (sensor_noise.position_error.size() != Nact ||
        sensor_noise.velocity_error.size() != Nact ||
        sensor_noise.acceleration_error.size() != Nact) {
        throw std::invalid_argument("Invalid sensor noise format");
    }
}

}; // namespace RAPTOR