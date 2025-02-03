#include "AccelerationLimits.h"

namespace RAPTOR {

AccelerationLimits::AccelerationLimits(std::shared_ptr<Trajectories>& trajPtr_input,
                                       const VecX& lowerLimits_input, 
                                       const VecX& upperLimits_input) : 
    trajPtr_(trajPtr_input),
    lowerLimits(lowerLimits_input), 
    upperLimits(upperLimits_input) {
    if (lowerLimits.size() != upperLimits.size()) {
        throw std::invalid_argument("lowerLimits and upperLimits must be the same size");
    }

    if (lowerLimits.size() != trajPtr_->Nact) {
        throw std::invalid_argument("lowerLimits and upperLimits must be the same size as the number of actuated joints");
    }

    initialize_memory(trajPtr_->N * trajPtr_->Nact, 
                      trajPtr_->varLength);
}

void AccelerationLimits::compute(const VecX& z, 
                                 bool compute_derivatives,
                                 bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    for (int i = 0; i < trajPtr_->N; i++) {
        g.segment(i * trajPtr_->Nact, trajPtr_->Nact) = trajPtr_->q_dd(i).head(trajPtr_->Nact);

        if (compute_derivatives) {
            pg_pz.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, trajPtr_->varLength) = trajPtr_->pq_dd_pz(i);
        }

        if (compute_hessian) {
            for (int j = 0; j < trajPtr_->Nact; j++) {
                pg_pz_pz(i * trajPtr_->Nact + j) = trajPtr_->pq_dd_pz_pz(j, i);
            }
        }
    }
}

void AccelerationLimits::compute_bounds() {
    for (int i = 0; i < trajPtr_->N; i++) {
        g_lb.segment(i * trajPtr_->Nact, trajPtr_->Nact) = lowerLimits;
        g_ub.segment(i * trajPtr_->Nact, trajPtr_->Nact) = upperLimits;
    }
}

void AccelerationLimits::print_violation_info() {
    for (int i = 0; i < trajPtr_->N; i++) {
        for (int j = 0; j < trajPtr_->Nact; j++) {
            if (g(i * trajPtr_->Nact + j) <= g_lb(i * trajPtr_->Nact + j)) {
                std::cout << "        AccelerationLimits.cpp: Joint " 
                          << j 
                          << " at time instance " 
                          << i 
                          << " is below lower limit: " 
                          << g(i * trajPtr_->Nact + j) 
                          << " < " 
                          << lowerLimits(j) 
                          << std::endl;
            }
            if (g(i * trajPtr_->Nact + j) >= g_ub(i * trajPtr_->Nact + j)) {
                std::cout << "        AccelerationLimits.cpp: Joint " 
                          << j 
                          << " at time instance "
                          << i 
                          << " is above upper limit: " 
                          << g(i * trajPtr_->Nact + j) 
                          << " > " 
                          << upperLimits(j) 
                          << std::endl;
            }
        }
    }
}

}; // namespace RAPTOR
