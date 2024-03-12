#include "JointLimits.h"

namespace IDTO {

JointLimits::JointLimits(std::shared_ptr<Trajectories>& trajPtr_input,
                         const VecX& lowerLimits_input, 
                         const VecX& upperLimits_input,
                         const bool wrapToPiOrNot_input) : 
    lowerLimits(lowerLimits_input), 
    upperLimits(upperLimits_input),
    wrapToPiOrNot(wrapToPiOrNot_input) {
    trajPtr_ = trajPtr_input;

    if (lowerLimits.size() != upperLimits.size()) {
        throw std::invalid_argument("lowerLimits and upperLimits must be the same size");
    }

    if (lowerLimits.size() != trajPtr_->Nact) {
        throw std::invalid_argument("lowerLimits and upperLimits must be the same size as the number of actuated joints");
    }

    m = trajPtr_->N * trajPtr_->Nact;
    
    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, trajPtr_->varLength);
}

void JointLimits::compute(const VecX& z, bool compute_derivatives) {
    if (is_computed(z, compute_derivatives)) {
        return;
    }

    if (compute_derivatives) {
        pg_pz.setZero();
    }

    trajPtr_->compute(z, compute_derivatives);

    for (int i = 0; i < trajPtr_->N; i++) {
        g.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, 1) = wrapToPiOrNot ? 
        							wrapToPi(trajPtr_->q(i).head(trajPtr_->Nact)) :
        							trajPtr_->q(i).head(trajPtr_->Nact);

        if (compute_derivatives) {
            pg_pz.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, trajPtr_->varLength) = trajPtr_->pq_pz(i);
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
            if (g(i * trajPtr_->Nact + j) <= g_lb(i * trajPtr_->Nact + j)) {
                std::cout << "        JointLimits.cpp: Actuator " << j 
                          << " at time instance " << i 
                          << " is violating the lower position limit" 
                          << std::endl;
            } 
            else if (g(i * trajPtr_->Nact + j) >= g_ub(i * trajPtr_->Nact + j)) {
                std::cout << "        JointLimits.cpp: Actuator " << j 
                          << " at time instance " << i 
                          << " is violating the upper position limit" 
                          << std::endl;
            }
        }
    }
}

}; // namespace IDTO
