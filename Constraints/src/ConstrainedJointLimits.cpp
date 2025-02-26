#include "ConstrainedJointLimits.h"

namespace RAPTOR {

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

void ConstrainedJointLimits::compute(const VecX& z, 
                                     bool compute_derivatives,
                                     bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    if (compute_hessian) {
        throw std::invalid_argument("ConstrainedJointLimits does not support hessian computation");
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    for (int i = 0; i < trajPtr_->N; i++) {
        VecX qfull(NB);
        dcPtr_->fill_independent_vector(qfull, trajPtr_->q(i), true);

        dcPtr_->setupJointPosition(qfull, compute_derivatives);
        g.segment(i * NB, NB) = Utils::wrapToPi(qfull);

        if (compute_derivatives) {
            // fill in independent joints derivatives directly
            for (int j = 0; j < dcPtr_->numIndependentJoints; j++) {
                int indenpendentJointIndex = dcPtr_->return_independent_joint_index(j);
                pg_pz.row(i * NB + indenpendentJointIndex) = trajPtr_->pq_pz(i).row(j);
            }

            // compute and fill in dependent joints derivatives
            MatX pq_dep_pz = dcPtr_->pq_dep_pq_indep * trajPtr_->pq_pz(i);
            for (int j = 0; j < dcPtr_->numDependentJoints; j++) {
                int denpendentJointIndex = dcPtr_->return_dependent_joint_index(j);
                pg_pz.row(i * NB + denpendentJointIndex) = pq_dep_pz.row(j);
            }
        }
    }
}

void ConstrainedJointLimits::compute_bounds() {
    for (int i = 0; i < trajPtr_->N; i++) {
        g_lb.segment(i * NB, NB) = lowerLimits;
        g_ub.segment(i * NB, NB) = upperLimits;
    }
}

void ConstrainedJointLimits::print_violation_info() {
    for (int i = 0; i < trajPtr_->N; i++) {
        for (int j = 0; j < NB; j++) {
            if (g(i * NB + j) <= lowerLimits(j)) {
                std::cout << "        ConstrainedJointLimits.cpp: Joint " 
                          << j 
                          << " at time instance " 
                          << i
                          << " is below lower limit: "
                          << g(i * NB + j) 
                          << " < " 
                          << lowerLimits(j) 
                          << std::endl;
            } 
            else if (g(i * NB + j) >= upperLimits(j)) {
                std::cout << "        ConstrainedJointLimits.cpp: Joint " 
                          << j 
                          << " at time instance " 
                          << i 
                          << " is above upper limit: " 
                          << g(i * NB + j) 
                          << " > " 
                          << upperLimits(j) 
                          << std::endl;
            }
        }
    }
}

}; // namespace RAPTOR
