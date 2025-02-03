#include "TrajectoryData.h"

namespace RAPTOR {

TrajectoryData::TrajectoryData(const std::string& filename_input,
                               const SensorNoiseInfo sensor_noise_input,
                               const TimeFormat time_format,
                               const int downsample_rate,
                               const bool add_sensor_noise) :
    sensor_noise(sensor_noise_input) {
    // parse file
    MatX traj_data = Utils::initializeEigenMatrixFromFile(filename_input);

    if ((traj_data.cols() - 1) % 3 == 0) {
        Nact = (traj_data.cols() - 1) / 3;
    }
    else if ((traj_data.cols() - 1) % 2 == 0) {
        Nact = (traj_data.cols() - 1) / 2;
    }
    else {
        std::cerr << "The file format needs to be:\n"
                  << "    [time, positions, velocities, accelerations/torques] or\n"
                  << "    [time, positions, velocities]" 
                  << std::endl;
        throw std::invalid_argument("Invalid trajectory file format!");
    }

    if (Nact <= 0) {
        throw std::invalid_argument("Invalid trajectory file format! 0 actuated joints!");
    }

    if (2 * Nact + 1 != traj_data.cols() &&
        3 * Nact + 1 != traj_data.cols()) {
        std::cerr << Nact << ' ' << traj_data.cols() << std::endl;
        throw std::invalid_argument("Invalid trajectory file format");
    }

    N = traj_data.rows();

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
        N = traj_data.rows();
    }

    tspan = traj_data.col(0);

    // check time format
    if (time_format == TimeFormat::Millisecond) {
        tspan *= 1e-3;
    }
    else if (time_format == TimeFormat::Microsecond) {
        tspan *= 1e-6;
    }
    else if (time_format == TimeFormat::Nanosecond) {
        tspan *= 1e-9;
    }

    for (int i = 1; i < tspan.size(); i++) {
        if (tspan(i - 1) > tspan(i)) {
            throw std::invalid_argument("Invalid time format: time is not increasing");
        }
    }

    T = tspan(N - 1) - tspan(0);

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

    // add sensor noise
    if (add_sensor_noise) {
        if (sensor_noise.position_error_type == SensorNoiseInfo::SensorNoiseType::Constant) {
            for (int i = 0; i < N; i++) {
                VecX noise = VecX::Random(Nact).cwiseProduct(sensor_noise.position_error);
                q(i) += noise;
            }
        }
        else if (sensor_noise.position_error_type == SensorNoiseInfo::SensorNoiseType::Ratio) {
            for (int i = 0; i < N; i++) {
                VecX noise = VecX::Random(Nact).cwiseProduct(q(i).cwiseProduct(sensor_noise.position_error));
                q(i) += q(i).cwiseProduct(sensor_noise.position_error);
            }
        }
        if (sensor_noise.velocity_error_type == SensorNoiseInfo::SensorNoiseType::Constant) {
            for (int i = 0; i < N; i++) {
                VecX noise = VecX::Random(Nact).cwiseProduct(sensor_noise.velocity_error);
                q_d(i) += noise;
            }
        }
        else if (sensor_noise.velocity_error_type == SensorNoiseInfo::SensorNoiseType::Ratio) {
            for (int i = 0; i < N; i++) {
                VecX noise = VecX::Random(Nact).cwiseProduct(q_d(i).cwiseProduct(sensor_noise.velocity_error));
                q_d(i) += q_d(i).cwiseProduct(sensor_noise.velocity_error);
            }
        }
        if (sensor_noise.acceleration_error_type == SensorNoiseInfo::SensorNoiseType::Constant) {
            for (int i = 0; i < N; i++) {
                VecX noise = VecX::Random(Nact).cwiseProduct(sensor_noise.acceleration_error);
                q_dd(i) += noise;
            }
        }
        else if (sensor_noise.acceleration_error_type == SensorNoiseInfo::SensorNoiseType::Ratio) {
            for (int i = 0; i < N; i++) {
                VecX noise = VecX::Random(Nact).cwiseProduct(q_dd(i).cwiseProduct(sensor_noise.acceleration_error));
                q_dd(i) += q_dd(i).cwiseProduct(sensor_noise.acceleration_error);
            }
        }
    }
}

TrajectoryData::TrajectoryData(const std::vector<std::string>& filenames_input,
                               const SensorNoiseInfo sensor_noise_input,
                               const TimeFormat time_format,
                               const int downsample_rate,
                               const bool add_sensor_noise) :
    sensor_noise(sensor_noise_input) {
    // parse file
    std::vector<MatX> traj_datas;
    N = 0;
    size_t NumCols = 0;
    for (const auto& filename : filenames_input) {
        try {
            traj_datas.push_back(Utils::initializeEigenMatrixFromFile(filename));
        }
        catch (const std::exception& e) {
            std::cerr << "Failed to load the file: " << filename << std::endl;
            std::cerr << "    " << e.what() << std::endl;
        }
        N += traj_datas.back().rows();
        if (NumCols == 0) {
            NumCols = traj_datas.back().cols();
        }
        else if (NumCols != traj_datas.back().cols()) {
            std::cerr << filename << std::endl;
            std::cerr << NumCols << ' ' << traj_datas.back().cols() << std::endl;
            throw std::invalid_argument("The trajectory files have different number of columns!");
        }
    }

    if (N == 0) {
        throw std::invalid_argument("No data loaded from the files!");
    }

    if ((NumCols - 1) % 3 == 0) {
        Nact = (NumCols - 1) / 3;
    }
    else if ((NumCols - 1) % 2 == 0) {
        Nact = (NumCols - 1) / 2;
    }
    else {
        std::cerr << "The file format needs to be:\n"
                  << "    [time, positions, velocities, accelerations/torques] or\n"
                  << "    [time, positions, velocities]" 
                  << std::endl;
        throw std::invalid_argument("Invalid trajectory file format!");
    }

    if (Nact <= 0) {
        throw std::invalid_argument("Invalid trajectory file format! 0 actuated joints!");
    }

    if (2 * Nact + 1 != NumCols &&
        3 * Nact + 1 != NumCols) {
        std::cerr << Nact << ' ' << NumCols << std::endl;
        throw std::invalid_argument("Invalid trajectory file format");
    }

    // check downsample format
    if (downsample_rate <= 0) {
        throw std::invalid_argument("Invalid downsample rate");
    }
    else if (downsample_rate > 1) {
        // int num_samples = traj_data.rows() / downsample_rate;
        // std::cout << "Performing downsample from " << traj_data.rows() << " to " << num_samples << std::endl;
        // MatX new_traj_data(num_samples, traj_data.cols());
        // Utils::uniformlySampleMatrixInRows(traj_data, new_traj_data, num_samples);
        // traj_data = new_traj_data; // this operation can not be reversed!
        // N = traj_data.rows();
    }
    
    std::vector<double> tspan_data;
    tspan_data.reserve(N);
    for (int i = 0; i < traj_datas.size(); i++) {
        for (int j = 0; j < traj_datas[i].rows(); j++) {
            tspan_data.push_back(traj_datas[i](j, 0));
        }
    }

    tspan = Eigen::Map<Eigen::VectorXd>(tspan_data.data(), tspan_data.size());

    // check time format
    if (time_format == TimeFormat::Millisecond) {
        tspan *= 1e-3;
    }
    else if (time_format == TimeFormat::Microsecond) {
        tspan *= 1e-6;
    }
    else if (time_format == TimeFormat::Nanosecond) {
        tspan *= 1e-9;
    }

    for (int i = 1; i < tspan.size(); i++) {
        if (tspan(i - 1) > tspan(i)) {
            throw std::invalid_argument("Invalid time format: time is not increasing!");
        }
    }

    T = tspan(N - 1) - tspan(0);

    // output information
    std::cout << "TrajectoryData: " << N << " data points loaded from a sequence of files" << std::endl;
    std::cout << "TrajectoryData: T = " << T << std::endl;
    std::cout << "TrajectoryData: Nact = " << Nact << std::endl;

    // no variable length parameters, we don't compute gradient here
    // trajectory data is usually large so save some memory here
    varLength = 0;
    initialize_memory();

    size_t index = 0;
    for (int i = 0; i < traj_datas.size(); i++) {
        for (int j = 0; j < traj_datas[i].rows(); j++) {
            q(index) = traj_datas[i].row(j).segment(1, Nact);
            q_d(index) = traj_datas[i].row(j).segment(1 + Nact, Nact);
            index++;
        }
    }

    // this could be acceleration estimation or applied torque from the sensor
    if (3 * Nact + 1 == NumCols) {
        std::cout << "TrajectoryData: q_dd or torque is loaded from the file" << std::endl;
        index = 0;
        for (int i = 0; i < traj_datas.size(); i++) {
            for (int j = 0; j < traj_datas[i].rows(); j++) {
                q_dd(index) = traj_datas[i].row(j).segment(1 + 2 * Nact, Nact);
                index++;
            }
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

    // add sensor noise
    if (add_sensor_noise) {
        if (sensor_noise.position_error_type == SensorNoiseInfo::SensorNoiseType::Constant) {
            for (int i = 0; i < N; i++) {
                VecX noise = VecX::Random(Nact).cwiseProduct(sensor_noise.position_error);
                q(i) += noise;
            }
        }
        else if (sensor_noise.position_error_type == SensorNoiseInfo::SensorNoiseType::Ratio) {
            for (int i = 0; i < N; i++) {
                VecX noise = VecX::Random(Nact).cwiseProduct(q(i).cwiseProduct(sensor_noise.position_error));
                q(i) += q(i).cwiseProduct(sensor_noise.position_error);
            }
        }
        if (sensor_noise.velocity_error_type == SensorNoiseInfo::SensorNoiseType::Constant) {
            for (int i = 0; i < N; i++) {
                VecX noise = VecX::Random(Nact).cwiseProduct(sensor_noise.velocity_error);
                q_d(i) += noise;
            }
        }
        else if (sensor_noise.velocity_error_type == SensorNoiseInfo::SensorNoiseType::Ratio) {
            for (int i = 0; i < N; i++) {
                VecX noise = VecX::Random(Nact).cwiseProduct(q_d(i).cwiseProduct(sensor_noise.velocity_error));
                q_d(i) += q_d(i).cwiseProduct(sensor_noise.velocity_error);
            }
        }
        if (sensor_noise.acceleration_error_type == SensorNoiseInfo::SensorNoiseType::Constant) {
            for (int i = 0; i < N; i++) {
                VecX noise = VecX::Random(Nact).cwiseProduct(sensor_noise.acceleration_error);
                q_dd(i) += noise;
            }
        }
        else if (sensor_noise.acceleration_error_type == SensorNoiseInfo::SensorNoiseType::Ratio) {
            for (int i = 0; i < N; i++) {
                VecX noise = VecX::Random(Nact).cwiseProduct(q_dd(i).cwiseProduct(sensor_noise.acceleration_error));
                q_dd(i) += q_dd(i).cwiseProduct(sensor_noise.acceleration_error);
            }
        }
    }
}

TrajectoryData::TrajectoryData(double T_input,
                               int N_input, 
                               int Nact_input,
                               const bool random_or_not,
                               const SensorNoiseInfo sensor_noise_input) :
    Trajectories(0, T_input, N_input, Nact_input, TimeDiscretization::Uniform),
    sensor_noise(sensor_noise_input) {
    if (random_or_not) {
        // randomly generate trajectory data
        for (int i = 0; i < N; i++) {
            q(i).setRandom();
            q_d(i).setRandom();
            q_dd(i).setRandom();
        }
    }
    else {
        // zero trajectory data
        // already did this in the constructor of Trajectories
        // for (int i = 0; i < N; i++) {
        //     q(i).setZero();
        //     q_d(i).setZero();
        //     q_dd(i).setZero();
        // }
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