#include "MinimizeInitialAcceleration.h"

namespace RAPTOR {

MinimizeInitialAcceleration::MinimizeInitialAcceleration(std::shared_ptr<Trajectories>& trajPtr_input) :
    trajPtr_(trajPtr_input) {
    initialize_memory(trajPtr_->varLength);
}

void MinimizeInitialAcceleration::compute(const VecX& z, 
                                          bool compute_derivatives,
                                          bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    // this will be called in idPtr_->compute
    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    const VecX& a0 = trajPtr_->q_dd(0);
    const double a0_squared = a0.dot(a0);
    const double a0_norm = std::sqrt(a0_squared);
    f = a0_norm;

    if (a0_norm > InitialAcceleration::SQUARE_ROOT_THRESHOLD) {
        if (compute_derivatives) {
            const MatX& pa0_pz = trajPtr_->pq_dd_pz(0);
            const VecX pa0_square_pz = pa0_pz.transpose() * a0;
            grad_f = pa0_square_pz / a0_norm;

            if (compute_hessian) {
                hess_f = pa0_pz.transpose() * pa0_pz / a0_norm;
                for (int j = 0; j < a0.size(); j++) {
                    const MatX& pa0_pz_pz = trajPtr_->pq_dd_pz_pz(j, 0);
                    hess_f += pa0_pz_pz * a0(j) / a0_norm;
                }
                hess_f -= pa0_square_pz * pa0_square_pz.transpose() / std::pow(a0_squared, 1.5);
            }
        }
    }
}

}; // namespace RAPTOR