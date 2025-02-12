#include "TrajectoryTerminalConstraints.h"

namespace RAPTOR {

TrajectoryTerminalConstraints::TrajectoryTerminalConstraints(std::shared_ptr<Trajectories>& trajPtr_input, 
                                                             const VecX desiredPoistion_input,
                                                             const VecX desiredVelocity_input,
                                                             const VecX desiredAcceleration_input) :
    trajPtr_(trajPtr_input),
    desiredPosition(desiredPoistion_input),
    desiredVelocity(desiredVelocity_input),
    desiredAcceleration(desiredAcceleration_input) {
    size_t m = 0;

    if (desiredPosition.size() == 0) {
        constrainTerminalPosition = false;
    }
    else if (desiredPosition.size() != trajPtr_->Nact) {
        throw std::invalid_argument("Error: desired position size is not equal to the number of actuators.");
    }
    else {
        constrainTerminalPosition = true;
        m += trajPtr_->Nact;
    }

    if (desiredVelocity.size() == 0) {
        constrainTerminalVelocity = false;
    }
    else if (desiredVelocity.size() != trajPtr_->Nact) {
        throw std::invalid_argument("Error: desired velocity size is not equal to the number of actuators.");
    }
    else {
        constrainTerminalVelocity = true;
        m += trajPtr_->Nact;
    }

    if (desiredAcceleration.size() == 0) {
        constrainTerminalAcceleration = false;
    }
    else if (desiredAcceleration.size() != trajPtr_->Nact) {
        throw std::invalid_argument("Error: desired acceleration size is not equal to the number of actuators.");
    }
    else {
        constrainTerminalAcceleration = true;
        m += trajPtr_->Nact;
    }

    initialize_memory(m, trajPtr_->varLength);
}

void TrajectoryTerminalConstraints::compute(const VecX& z, 
                                            bool compute_derivatives,
                                            bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    size_t g_index = 0;
    const size_t last_time_instane = trajPtr_->N - 1;

    if (constrainTerminalPosition) {
        g.segment(g_index, trajPtr_->Nact) = trajPtr_->q(last_time_instane);

        if (compute_derivatives) {
            pg_pz.block(g_index, 0, trajPtr_->Nact, trajPtr_->varLength) = trajPtr_->pq_pz(last_time_instane);
        }

        if (compute_hessian) {
            for (int j = 0; j < trajPtr_->Nact; j++) {
                pg_pz_pz(g_index + j) = trajPtr_->pq_pz_pz(j, last_time_instane);
            }
        }

        g_index += trajPtr_->Nact;
    }

    if (constrainTerminalVelocity) {
        g.segment(g_index, trajPtr_->Nact) = trajPtr_->q_d(last_time_instane);

        if (compute_derivatives) {
            pg_pz.block(g_index, 0, trajPtr_->Nact, trajPtr_->varLength) = trajPtr_->pq_d_pz(last_time_instane);
        }

        if (compute_hessian) {
            for (int j = 0; j < trajPtr_->Nact; j++) {
                pg_pz_pz(g_index + j) = trajPtr_->pq_d_pz_pz(j, last_time_instane);
            }
        }

        g_index += trajPtr_->Nact;
    }

    if (constrainTerminalAcceleration) {
        g.segment(g_index, trajPtr_->Nact) = trajPtr_->q_dd(last_time_instane);

        if (compute_derivatives) {
            pg_pz.block(g_index, 0, trajPtr_->Nact, trajPtr_->varLength) = trajPtr_->pq_dd_pz(last_time_instane);
        }

        if (compute_hessian) {
            for (int j = 0; j < trajPtr_->Nact; j++) {
                pg_pz_pz(g_index + j) = trajPtr_->pq_dd_pz_pz(j, last_time_instane);
            }
        }

        g_index += trajPtr_->Nact;
    }
}

void TrajectoryTerminalConstraints::compute_bounds() {
    size_t g_index = 0;

    if (constrainTerminalPosition) {
        g_lb.segment(g_index, trajPtr_->Nact) = desiredPosition;
        g_ub.segment(g_index, trajPtr_->Nact) = desiredPosition;
        g_index += trajPtr_->Nact;
    }

    if (constrainTerminalVelocity) {
        g_lb.segment(g_index, trajPtr_->Nact) = desiredVelocity;
        g_ub.segment(g_index, trajPtr_->Nact) = desiredVelocity;
        g_index += trajPtr_->Nact;
    }

    if (constrainTerminalAcceleration) {
        g_lb.segment(g_index, trajPtr_->Nact) = desiredAcceleration;
        g_ub.segment(g_index, trajPtr_->Nact) = desiredAcceleration;
        g_index += trajPtr_->Nact;
    }
}

void TrajectoryTerminalConstraints::print_violation_info() {
    size_t g_index = 0;

    if (constrainTerminalPosition) {
        for (int j = 0; j < trajPtr_->Nact; j++) {
            if (g(g_index + j) <= g_lb(g_index + j) - 1e-4 ||
                g(g_index + j) >= g_ub(g_index + j) + 1e-4) {
                std::cout << "        TrajectoryTerminalConstraints.cpp: Terminal position constraint violation: " 
                          << g(g_index + j) 
                          << " != " 
                          << g_lb(g_index + j) 
                          << std::endl;

            }
        }

        g_index += trajPtr_->Nact;
    }

    if (constrainTerminalVelocity) {
        for (int j = 0; j < trajPtr_->Nact; j++) {
            if (g(g_index + j) <= g_lb(g_index + j) - 1e-4 ||
                g(g_index + j) >= g_ub(g_index + j) + 1e-4) {
                std::cout << "        TrajectoryTerminalConstraints.cpp: Terminal velocity constraint violation: " 
                          << g(g_index + j) 
                          << " != " 
                          << g_lb(g_index + j) 
                          << std::endl;
            }
        }

        g_index += trajPtr_->Nact;
    }

    if (constrainTerminalAcceleration) {
        for (int j = 0; j < trajPtr_->Nact; j++) {
            if (g(g_index + j) <= g_lb(g_index + j) - 1e-4 ||
                g(g_index + j) >= g_ub(g_index + j) + 1e-4) {
                std::cout << "        TrajectoryTerminalConstraints.cpp: Terminal acceleration constraint violation: " 
                          << g(g_index + j) 
                          << " != " 
                          << g_lb(g_index + j) 
                          << std::endl;
            }
        }

        g_index += trajPtr_->Nact;
    }
}

}; // namespace RAPTOR