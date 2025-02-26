#include "MinimizePower.h"

namespace RAPTOR {

MinimizePower::MinimizePower(std::shared_ptr<Trajectories>& trajPtr_input, 
                             std::shared_ptr<InverseDynamics>& idPtr_input) :
    trajPtr_(trajPtr_input),
    idPtr_(idPtr_input) {
    initialize_memory(trajPtr_->varLength);
}

void MinimizePower::compute(const VecX& z, 
                            bool compute_derivatives,
                            bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    // this will be called in idPtr_->compute
    // trajPtr_->compute(z, compute_derivatives, compute_hessian);

    idPtr_->compute(z, compute_derivatives, compute_hessian);

    f = 0;
    if (compute_derivatives) {
        grad_f.setZero();
        if (compute_hessian) {
            hess_f.setZero();
        }
    }

    for (int i = 0; i < trajPtr_->N; i++) {
        const VecX& tau = idPtr_->tau(i);
        const VecX& v = trajPtr_->q_d(i);
        const VecX power = tau.cwiseProduct(v);
        const double power_squared = power.dot(power);

        if (power_squared > Power::SQUARE_ROOT_THRESHOLD) {
            const double power_norm = std::sqrt(power_squared);
            f += power_norm;

            if (compute_derivatives) {
                const MatX& ptau_pz = idPtr_->ptau_pz(i);
                const MatX& pv_pz = trajPtr_->pq_d_pz(i);
                const MatX ppower_pz = tau.asDiagonal() * pv_pz + v.asDiagonal() * ptau_pz;
                const VecX ppower_square_pz = ppower_pz.transpose() * power;
                grad_f += ppower_square_pz / power_norm;

                if (compute_hessian) {
                    hess_f += ppower_pz.transpose() * ppower_pz / power_norm;
                    for (int j = 0; j < power.size(); j++) {
                        const MatX ppower_pz_pz = v(j) * idPtr_->ptau_pz_pz(j, i) +           // v * ∂²τ/∂z∂z
                                                  ptau_pz.row(j).transpose() * pv_pz.row(j) + // ∂τ/∂z * ∂v/∂z^T
                                                  pv_pz.row(j).transpose() * ptau_pz.row(j) + // ∂v/∂z * ∂τ/∂z^T
                                                  tau(j) * trajPtr_->pq_d_pz_pz(j, i);        // τ * ∂²v/∂z∂z
                        hess_f += ppower_pz_pz * power(j) / power_norm;
                    }
                    hess_f -= ppower_square_pz * ppower_square_pz.transpose() / std::pow(power_squared, 1.5);
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