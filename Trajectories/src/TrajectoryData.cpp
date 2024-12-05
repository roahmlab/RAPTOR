#include "TrajectoryData.h"

namespace RAPTOR {

TrajectoryData::TrajectoryData(const std::string& filename_input) {
    // parse file
    const MatX traj_data = Utils::initializeEigenMatrixFromFile(filename_input);

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
}

TrajectoryData::TrajectoryData(double T_input,
                               int N_input, 
                               int Nact_input) :
    Trajectories(0, T_input, N_input, Nact_input, TimeDiscretization::Uniform) {
    // randomly generate trajectory data
    for (int i = 0; i < N; i++) {
        q(i).setRandom();
        q_d(i).setRandom();
        q_dd(i).setRandom();
    }
}

}; // namespace RAPTOR