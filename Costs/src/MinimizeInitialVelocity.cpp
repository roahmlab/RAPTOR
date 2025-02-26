#include "MinimizeInitialVelocity.h"

namespace RAPTOR {

MinimizeInitialVelocity::MinimizeInitialVelocity(std::shared_ptr<Trajectories>& trajPtr_input) :
    trajPtr_(trajPtr_input) {
    initialize_memory(trajPtr_->varLength);
}

void MinimizeInitialVelocity::compute(const VecX& z, 
                                      bool compute_derivatives,
                                      bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    // this will be called in idPtr_->compute
    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    const VecX& v0 = trajPtr_->q_d(0);
    const double v0_squared = v0.dot(v0);
    const double v0_norm = std::sqrt(v0_squared);
    f = v0_norm;

    if (v0_norm > InitialVelocity::SQUARE_ROOT_THRESHOLD) {
        if (compute_derivatives) {
            const MatX& pv0_pz = trajPtr_->pq_d_pz(0);
            const VecX pv0_square_pz = pv0_pz.transpose() * v0;
            grad_f = pv0_square_pz / v0_norm;

            if (compute_hessian) {
                hess_f = pv0_pz.transpose() * pv0_pz / v0_norm;
                for (int j = 0; j < v0.size(); j++) {
                    const MatX& pv0_pz_pz = trajPtr_->pq_d_pz_pz(j, 0);
                    hess_f += pv0_pz_pz * v0(j) / v0_norm;
                }
                hess_f -= pv0_square_pz * pv0_square_pz.transpose() / std::pow(v0_squared, 1.5);
            }
        }
    }
}

}; // namespace RAPTOR