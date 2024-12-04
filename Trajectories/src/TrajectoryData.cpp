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
   
    if (3 * Nact + 1 != traj_data.cols()) {
        throw std::invalid_argument("Invalid trajectory file format");
    }

    tspan = traj_data.col(0)/1e9;
    T = tspan(N - 1);

    // output information
    std::cout << "TrajectoryData: " << N << " data points loaded from " << filename_input << std::endl;
    std::cout << "TrajectoryData: T = " << T << std::endl;
    std::cout << "TrajectoryData: Nact = " << Nact << std::endl;

    // load trajectories
    varLength = 0; // no variable length parameters, we don't compute gradient here
    initialize_memory();

    for (int i = 0; i < N; i++) {
        q(i) = traj_data.row(i).segment(1, Nact);
        q_d(i) = traj_data.row(i).segment(1 + Nact, Nact);
        q_dd(i).setZero();
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