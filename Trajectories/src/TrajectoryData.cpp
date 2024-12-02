#include "TrajectoryData.h"

namespace RAPTOR {

TrajectoryData::TrajectoryData(const std::string& filename_input) {
    // parse file
    const MatX traj_data = Utils::initializeEigenMatrixFromFile(filename_input);

    N = traj_data.rows();
    Nact = (traj_data.cols() - 1) / 2;

    if (Nact <= 0) {
        throw std::invalid_argument("0 actuated joints");
    }

    if (2 * Nact + 1 != traj_data.cols()) {
        throw std::invalid_argument("Invalid trajectory file format");
    }

    tspan = traj_data.col(0);
    T = tspan(N - 1);

    // output information
    std::cout << "TrajectoryData: " << N << " data points loaded from " << filename_input << std::endl;
    std::cout << "TrajectoryData: T = " << T << std::endl;
    std::cout << "TrajectoryData: Nact = " << Nact << std::endl;

    // load trajectories
    varLength = 3 * Nact; // decision variables for q, q_d, q_dd are just themselves
    initialize_memory();

    for (int i = 0; i < N; i++) {
        q(i) = traj_data.row(i).segment(1, Nact);
        q_d(i) = traj_data.row(i).segment(1 + Nact, Nact);
        q_dd(i).setZero();

        pq_pz(i).setZero();
        pq_d_pz(i).setZero();
        pq_dd_pz(i).setZero();

        pq_pz(i).middleCols(0, Nact).setIdentity();
        pq_d_pz(i).middleCols(Nact, Nact).setIdentity();
        pq_dd_pz(i).middleCols(2 * Nact, Nact).setIdentity();
    }
}

TrajectoryData::TrajectoryData(double T_input,
                               int N_input, 
                               int Nact_input) :
    Trajectories(3 * Nact_input, T_input, N_input, Nact_input, TimeDiscretization::Uniform) {
    // randomly generate trajectory data
    for (int i = 0; i < N; i++) {
        q(i).setRandom();
        q_d(i).setRandom();
        q_dd(i).setRandom();

        pq_pz(i).setZero();
        pq_d_pz(i).setZero();
        pq_dd_pz(i).setZero();

        pq_pz(i).middleCols(0, Nact).setIdentity();
        pq_d_pz(i).middleCols(Nact, Nact).setIdentity();
        pq_dd_pz(i).middleCols(2 * Nact, Nact).setIdentity();
    }
}

}; // namespace RAPTOR