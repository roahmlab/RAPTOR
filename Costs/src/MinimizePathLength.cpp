#include "MinimizePathLength.h"

namespace RAPTOR {

MinimizePathLength::MinimizePathLength(std::shared_ptr<Trajectories>& trajPtr_input) :
    trajPtr_(trajPtr_input) {
    initialize_memory(trajPtr_->varLength);
}

void MinimizePathLength::compute(const VecX& z, 
                                 bool compute_derivatives,
                                 bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    // this will be called in idPtr_->compute
    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    f = 0;
    if (compute_derivatives) {
        grad_f.setZero();
        if (compute_hessian) {
            hess_f.setZero();
        }
    }

    for (int i = 0; i < trajPtr_->N - 1; i++) {
        const VecX& q1 = trajPtr_->q(i);
        const VecX& q2 = trajPtr_->q(i + 1);
        const VecX q_diff = q2 - q1;

        f += 0.5 * q_diff.dot(q_diff);

        if (compute_derivatives) {
            const MatX& pq1_pz = trajPtr_->pq_pz(i);
            const MatX& pq2_pz = trajPtr_->pq_pz(i + 1);
            const MatX pq_diff_pz = pq2_pz - pq1_pz;

            grad_f += q_diff.transpose() * pq_diff_pz;

            if (compute_hessian) {
                hess_f += pq_diff_pz.transpose() * pq_diff_pz;

                for (int j = 0; j < q_diff.size(); j++) {
                    const MatX& pq1_pz_pz = trajPtr_->pq_pz_pz(j, i);
                    const MatX& pq2_pz_pz = trajPtr_->pq_pz_pz(j, i + 1);
                    const MatX pq_diff_pz_pz = pq2_pz_pz - pq1_pz_pz;
                    hess_f += q_diff(j) * pq_diff_pz_pz;
                }
            }
        }
    }

    f /= trajPtr_->N;
    if (compute_derivatives) {
        grad_f /= trajPtr_->N;
        if (compute_hessian) {
            hess_f /= trajPtr_->N;
        }
    }
}

}; // namespace RAPTOR