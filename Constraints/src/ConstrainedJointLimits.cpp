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
    if (compute_derivatives) {
        pg_pz.setZero();
    }

    trajPtr_->compute(z, compute_derivatives);

    for (int i = 0; i < trajPtr_->N; i++) {
        VecX qfull(NB);
        dcPtr_->fill_independent_vector(qfull, trajPtr_->q(i));
        dcPtr_->setupJointPosition(qfull, compute_derivatives);
        g.block(i * NB, 0, NB, 1) = qfull;

        if (compute_derivatives) {
            // fill in independent joints derivatives directly
            for (int j = 0; j < dcPtr_->numIndependentJoints; j++) {
                int indenpendentJointIndex = dcPtr_->return_independent_joint_index(j);
                pg_pz.row(i * NB + indenpendentJointIndex) = trajPtr_->pq_pz(i).row(j);
            }

            // quick sanity check
            assert(dcPtr_->pq_unact_pq_act.cols() == dcPtr_->numIndependentJoints);
            assert(dcPtr_->pq_unact_pq_act.rows() == dcPtr_->numDependentJoints);

            // compute and fill in dependent joints derivatives
            SpaMatX pq_unact_pz = dcPtr_->pq_unact_pq_act * trajPtr_->pq_pz(i);
            for (int j = 0; j < dcPtr_->numDependentJoints; j++) {
                int denpendentJointIndex = dcPtr_->return_dependent_joint_index(j);
                pg_pz.row(i * NB + denpendentJointIndex) = pq_unact_pz.row(j);
            }
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
