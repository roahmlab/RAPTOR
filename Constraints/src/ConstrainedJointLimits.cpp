#include "ConstrainedJointLimits.h"

namespace IDTO {

ConstrainedJointLimits::ConstrainedJointLimits(std::shared_ptr<Trajectories>& trajPtr_input, 
                                               std::shared_ptr<DynamicsConstraints>& dcPtr_input, 
                                               const VecX& lowerLimits_input, 
                                               const VecX& upperLimits_input) {
    trajPtr_ = trajPtr_input;
    dcPtr_ = dcPtr_input;
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
        dcPtr_->fill_independent_vector(qfull, trajPtr_->q(i), true);

        dcPtr_->setupJointPosition(qfull, compute_derivatives);
        g.block(i * NB, 0, NB, 1) = qfull;

        if (compute_derivatives) {
            // fill in independent joints derivatives directly
            for (int j = 0; j < dcPtr_->numIndependentJoints; j++) {
                int indenpendentJointIndex = dcPtr_->return_independent_joint_index(j);
                pg_pz.row(i * NB + indenpendentJointIndex) = trajPtr_->pq_pz(i).row(j);
            }

            // quick sanity check
            if (dcPtr_->pq_unact_pq_act.cols() != dcPtr_->numIndependentJoints) {
                throw std::runtime_error("pq_unact_pq_act must have the same number of columns as the number of independent joints");
            }
            if (dcPtr_->pq_unact_pq_act.rows() != dcPtr_->numDependentJoints) {
                throw std::runtime_error("pq_unact_pq_act must have the same number of rows as the number of dependent joints");
            }

            // compute and fill in dependent joints derivatives
            MatX pq_unact_pz = dcPtr_->pq_unact_pq_act * trajPtr_->pq_pz(i);
            for (int j = 0; j < dcPtr_->numDependentJoints; j++) {
                int denpendentJointIndex = dcPtr_->return_dependent_joint_index(j);
                pg_pz.row(i * NB + denpendentJointIndex) = pq_unact_pz.row(j);
            }
        }
    }
}

void ConstrainedJointLimits::compute_bounds() {
    for (int i = 0; i < trajPtr_->N; i++) {
        g_lb.block(i * NB, 0, NB, 1) = lowerLimits;
        g_ub.block(i * NB, 0, NB, 1) = upperLimits;
    }
}

}; // namespace IDTO
