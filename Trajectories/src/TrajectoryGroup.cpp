#include "TrajectoryGroup.h"

namespace IDTO {

void TrajectoryGroup::add_trajectory(const std::string& name,    
                                     std::shared_ptr<Trajectories> trajectory) {
    trajectories[name] = trajectory;
    gather_trajectories_information();
}

void TrajectoryGroup::gather_trajectories_information(const bool print_info) {
    T = 0;
    N = 0;
    Nact = trajectories.begin()->second->Nact;
    varLength = 0;

    if (print_info) std::cout << "Trajectory group information: " << std::endl;
    
    size_t index = 0;
    for (const auto& it : trajectories) {
        T += it.second->T;

        trajectory_locations[it.first] = std::make_pair(N, N + it.second->N);
        N += it.second->N;

        if (Nact != it.second->Nact) {
            throw std::invalid_argument("The number of actuated joints in the trajectories are not the same.");
        }

        variable_locations[it.first] = std::make_pair(varLength, varLength + it.second->varLength);
        varLength += it.second->varLength;

        if (print_info) {
            std::cout << "Trajectory " << index << ": " << it.first << std::endl;
            std::cout << "    trajectory location: [" 
                      << trajectory_locations[it.first].first << ", " 
                      << trajectory_locations[it.first].second << "]" << std::endl;
            std::cout << "    variable location: ["
                      << variable_locations[it.first].first << ", "
                      << variable_locations[it.first].second << "]" << std::endl;
        }

        index++;
    }

    q.resize(1, N);
    q_d.resize(1, N);
    q_dd.resize(1, N);

    pq_pz.resize(1, N);
    pq_d_pz.resize(1, N);
    pq_dd_pz.resize(1, N);

    pq_pz_pz.resize(Nact, N);
    pq_d_pz_pz.resize(Nact, N);
    pq_dd_pz_pz.resize(Nact, N);

    for (int i = 0; i < N; i++) {
        q(i) = VecX::Zero(Nact);
        q_d(i) = VecX::Zero(Nact);
        q_dd(i) = VecX::Zero(Nact);

        pq_pz(i) = MatX::Zero(Nact, varLength);
        pq_d_pz(i) = MatX::Zero(Nact, varLength);
        pq_dd_pz(i) = MatX::Zero(Nact, varLength);

        for (int j = 0; j < Nact; j++) {
            pq_pz_pz(j, i) = MatX::Zero(varLength, varLength);
            pq_d_pz_pz(j, i) = MatX::Zero(varLength, varLength);
            pq_dd_pz_pz(j, i) = MatX::Zero(varLength, varLength);
        }
    }
}

void TrajectoryGroup::compute(const VecX& z, 
                              bool compute_derivatives,
                              bool compute_hessian) {
    for (const auto& it : trajectories) {
        const size_t trajectory_offset = trajectory_locations[it.first].first;
        const size_t variable_offset = variable_locations[it.first].first;

        const VecX& z_segment = z.segment(variable_offset, 
                                          variable_locations[it.first].second - variable_offset);
        it.second->compute(z_segment, compute_derivatives, compute_hessian);

        // copy the computed values of each trajectory to the local variables in TrajectoryGroup
        for (int i = 0; i < it.second->N; i++) {
            q(trajectory_offset + i) = it.second->q(i);
            q_d(trajectory_offset + i) = it.second->q_d(i);
            q_dd(trajectory_offset + i) = it.second->q_dd(i);

            if (compute_derivatives) {
                pq_pz(trajectory_offset + i)
                    .block(0, variable_offset, Nact, it.second->varLength) 
                        = it.second->pq_pz(i);
                pq_d_pz(trajectory_offset + i)
                    .block(0, variable_offset, Nact, it.second->varLength) 
                        = it.second->pq_d_pz(i);
                pq_dd_pz(trajectory_offset + i)
                    .block(0, variable_offset, Nact, it.second->varLength) 
                        = it.second->pq_dd_pz(i);
            }

            if (compute_hessian) {
                for (int j = 0; j < Nact; j++) {
                    pq_pz_pz(j, trajectory_offset + i)
                        .block(variable_offset, variable_offset, it.second->varLength, it.second->varLength) 
                            = it.second->pq_pz_pz(j, i);
                    pq_d_pz_pz(j, trajectory_offset + i)
                        .block(variable_offset, variable_offset, it.second->varLength, it.second->varLength) 
                            = it.second->pq_d_pz_pz(j, i);
                    pq_dd_pz_pz(j, trajectory_offset + i)
                        .block(variable_offset, variable_offset, it.second->varLength, it.second->varLength) 
                            = it.second->pq_dd_pz_pz(j, i);
                }
            }
        }
    }
}

}; // namespace IDTO