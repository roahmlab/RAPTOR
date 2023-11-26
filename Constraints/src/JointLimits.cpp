#include "JointLimits.h"

namespace IDTO {

JointLimits::JointLimits(std::unique_ptr<Trajectories> trajPtr_input,
                         const VecX& lowerLimits_input, 
                         const VecX& upperLimits_input) : 
    lowerLimits(lowerLimits_input), 
    upperLimits(upperLimits_input) {
    trajPtr_ = std::move(trajPtr_input);

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
    if (compute_derivatives) {
        pg_pz.setZero();
    }

    trajPtr_->compute(z, compute_derivatives);

    for (int i = 0; i < trajPtr_->N; i++) {
        g.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, 1) = trajPtr_->q(i);

        if (compute_derivatives) {
            // block operation is not supported for sparse matrices
            // pg_pz.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, trajPtr_->varLength) = trajPtr_->pq_pz(i);
            for (int j = 0; j < trajPtr_->Nact; j++) {
                pg_pz.row(i * trajPtr_->Nact + j) = trajPtr_->pq_pz(i).row(j);
            }
        }
    }
}

void JointLimits::compute_lb() {
    for (int i = 0; i < trajPtr_->N; i++) {
        g_lb.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, 1) = lowerLimits;
    }
}

void JointLimits::compute_ub() {
    for (int i = 0; i < trajPtr_->N; i++) {
        g_ub.block(i * trajPtr_->Nact, 0, trajPtr_->Nact, 1) = upperLimits;
    }
}

}; // namespace IDTO
