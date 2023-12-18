#include "TorqueLimits.h"

namespace IDTO {

TorqueLimits::TorqueLimits(std::shared_ptr<Trajectories>& trajPtr_input, 
                           std::shared_ptr<InverseDynamics> idPtr_input,
                           const VecX& lowerLimits_input, 
                           const VecX& upperLimits_input) :
    trajPtr_(trajPtr_input),
    idPtr_(idPtr_input) {
    lowerLimits = lowerLimits_input;
    upperLimits = upperLimits_input;
    
    m = trajPtr_->N * trajPtr_->Nact;

    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, trajPtr_->varLength);
}

void TorqueLimits::compute(const VecX& z, bool compute_derivatives) {
    if (is_computed(z, compute_derivatives)) {
        return;
    }

    if (compute_derivatives) {
        pg_pz.setZero();
    }

    trajPtr_->compute(z, compute_derivatives);

    idPtr_->compute(z, compute_derivatives);

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

}; // namespace IDTO