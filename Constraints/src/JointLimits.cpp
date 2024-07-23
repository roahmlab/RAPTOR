#include "JointLimits.h"

namespace RAPTOR {

JointLimits::JointLimits(std::shared_ptr<Trajectories>& trajPtr_input,
                         const VecX& lowerLimits_input, 
                         const VecX& upperLimits_input) : 
    lowerLimits(lowerLimits_input), 
    upperLimits(upperLimits_input) {
    trajPtr_ = trajPtr_input;

    if (lowerLimits.size() != upperLimits.size()) {
        throw std::invalid_argument("lowerLimits and upperLimits must be the same size");
    }

    if (lowerLimits.size() != trajPtr_->Nact) {
        throw std::invalid_argument("lowerLimits and upperLimits must be the same size as the number of actuated joints");
    }

    initialize_memory(trajPtr_->N * trajPtr_->Nact, 
                      trajPtr_->varLength);
}

void JointLimits::compute(const VecX& z, 
                          bool compute_derivatives,
                          bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    for (int i = 0; i < trajPtr_->N; i++) {
        // g.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, 1) = Utils::wrapToPi(trajPtr_->q(i).head(trajPtr_->Nact));
        g.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, 1) = trajPtr_->q(i).head(trajPtr_->Nact);

        if (compute_derivatives) {
            pg_pz.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, trajPtr_->varLength) = trajPtr_->pq_pz(i);
        }

        if (compute_hessian) {
            for (int j = 0; j < trajPtr_->Nact; j++) {
                pg_pz_pz(i * trajPtr_->Nact + j) = trajPtr_->pq_pz_pz(j, i);
            }
        }
    }
}

void JointLimits::compute_bounds() {
    for (int i = 0; i < trajPtr_->N; i++) {
        g_lb.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, 1) = lowerLimits;
        g_ub.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, 1) = upperLimits;
    }
}

void JointLimits::print_violation_info() {
    for (int i = 0; i < trajPtr_->N; i++) {
        for (int j = 0; j < trajPtr_->Nact; j++) {
            if (g(i * trajPtr_->Nact + j) <= lowerLimits(j)) {
                std::cout << "        JointLimits.cpp: Joint " 
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
                std::cout << "        JointLimits.cpp: Joint " 
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
