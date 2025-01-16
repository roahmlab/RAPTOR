#include "MinimizeJerk.h"

namespace RAPTOR {

MinimizeJerk::MinimizeJerk(std::shared_ptr<Trajectories>& trajPtr_input) :
    trajPtr_(trajPtr_input) {
    initialize_memory(trajPtr_->varLength);
}

void MinimizeJerk::compute(const VecX& z, 
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
        const VecX& q_dd1 = trajPtr_->q_dd(i);
        const VecX& q_dd2 = trajPtr_->q_dd(i + 1);
        const VecX q_dd_diff = q_dd2 - q_dd1;

        f += 0.5 * q_dd_diff.dot(q_dd_diff);

        if (compute_derivatives) {
            const MatX& pq_dd1_pz = trajPtr_->pq_dd_pz(i);
            const MatX& pq_dd2_pz = trajPtr_->pq_dd_pz(i + 1);
            const MatX pq_dd_diff_pz = pq_dd2_pz - pq_dd1_pz;

            grad_f += q_dd_diff.transpose() * pq_dd_diff_pz;

            if (compute_hessian) {
                hess_f += pq_dd_diff_pz.transpose() * pq_dd_diff_pz;

                for (int j = 0; j < q_dd_diff.size(); j++) {
                    const MatX& pq_dd1_pz_pz = trajPtr_->pq_dd_pz_pz(j, i);
                    const MatX& pq_dd2_pz_pz = trajPtr_->pq_dd_pz_pz(j, i + 1);
                    const MatX pq_dd_diff_pz_pz = pq_dd2_pz_pz - pq_dd1_pz_pz;
                    hess_f += q_dd_diff(j) * pq_dd_diff_pz_pz;
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