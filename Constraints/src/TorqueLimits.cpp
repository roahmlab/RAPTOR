#include "TorqueLimits.h"

namespace RAPTOR {

TorqueLimits::TorqueLimits(std::shared_ptr<Trajectories>& trajPtr_input, 
                           std::shared_ptr<InverseDynamics> idPtr_input,
                           const VecX& lowerLimits_input, 
                           const VecX& upperLimits_input) :
    trajPtr_(trajPtr_input),
    idPtr_(idPtr_input) {
    lowerLimits = lowerLimits_input;
    upperLimits = upperLimits_input;
    
    initialize_memory(trajPtr_->N * trajPtr_->Nact, 
                      trajPtr_->varLength,
                      false);
}

void TorqueLimits::compute(const VecX& z, 
                           bool compute_derivatives,
                           bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    if (compute_hessian) {
        throw std::invalid_argument("TorqueLimits: Hessian not implemented yet");
    }

    idPtr_->compute(z, compute_derivatives, compute_hessian);

    for (int i = 0; i < trajPtr_->N; i++) {
        g.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, 1) = idPtr_->tau(i);

        if (compute_derivatives) {
            pg_pz.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, trajPtr_->varLength) = idPtr_->ptau_pz(i);
        }
    }
}

void TorqueLimits::compute_bounds() {
    for (int i = 0; i < trajPtr_->N; i++) {
        g_lb.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, 1) = lowerLimits;
        g_ub.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, 1) = upperLimits;
    }
}

void TorqueLimits::print_violation_info() {
    for (int i = 0; i < trajPtr_->N; i++) {
        for (int j = 0; j < trajPtr_->Nact; j++) {
            if (g(i * trajPtr_->Nact + j) <= lowerLimits(j)) {
                std::cout << "        TorqueLimits.cpp: Joint " 
                          << j 
                          << " at time instance " 
                          << i
                          << " is below lower limit: "
                          << g(i * trajPtr_->Nact + j) 
                          << " < " 
                          << lowerLimits(j) 
                          << std::endl;
            } 
            else if (g(i * trajPtr_->Nact + j) >= upperLimits(j)) {
                std::cout << "        TorqueLimits.cpp: Joint " 
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