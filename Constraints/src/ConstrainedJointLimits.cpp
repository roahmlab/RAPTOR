#include "ConstrainedJointLimits.h"

namespace IDTO {

ConstrainedJointLimits::ConstrainedJointLimits(std::unique_ptr<Trajectories> trajPtr_input, 
                                               std::unique_ptr<DynamicsConstraints> dcPtr_input, 
                                               const VecX& lowerLimits_input, 
                                               const VecX& upperLimits_input) {
    trajPtr_ = std::move(trajPtr_input);
    dcPtr_ = std::move(dcPtr_input);
    lowerLimits = lowerLimits_input;
    upperLimits = upperLimits_input;

    if (lowerLimits.size() != upperLimits.size()) {
        throw std::invalid_argument("lowerLimits and upperLimits must be the same size");
    }

    if (trajPtr_->Nact != dcPtr_->numIndependentJoints) {
        throw std::invalid_argument("Trajectory and DynamicsConstraints must have the same number of actuated joints");
    }

    NB = dcPtr_->numIndependentJoints + dcPtr_->numDependentJoints;

    if (lowerLimits.size() != NB) {
        throw std::invalid_argument("lowerLimits and upperLimits must be the same size as the number of actuated joints");
    }

    m = trajPtr_->N * NB;
    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, trajPtr_->varLength);
}

void ConstrainedJointLimits::compute(const VecX& z, bool compute_derivatives) {
    trajPtr_->compute(z, compute_derivatives);

    for (int i = 0; i < trajPtr_->N; i++) {
        g.block(i * NB, 0, NB, 1) = trajPtr_->q(i);
    }

    if (compute_derivatives) {
        pg_pz.setZero();

        for (int i = 0; i < trajPtr_->N; i++) {
            // pg_pz.block(i * NB, 0, NB, varLength) = trajPtr_->pq_pz(i);
        }
    }
}

void ConstrainedJointLimits::compute_lb() {
    for (int i = 0; i < trajPtr_->N; i++) {
        g_lb.block(i * NB, 0, NB, 1) = lowerLimits;
    }
}

void ConstrainedJointLimits::compute_ub() {
    for (int i = 0; i < trajPtr_->N; i++) {
        g_ub.block(i * NB, 0, NB, 1) = upperLimits;
    }
}

}; // namespace IDTO
