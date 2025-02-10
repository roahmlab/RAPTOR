#include "MinimizeTorque.h"

namespace RAPTOR {

MinimizeTorque::MinimizeTorque(std::shared_ptr<Trajectories>& trajPtr_input, 
                               std::shared_ptr<InverseDynamics>& idPtr_input) :
    trajPtr_(trajPtr_input),
    idPtr_(idPtr_input) {
    initialize_memory(trajPtr_->varLength);
}

void MinimizeTorque::compute(const VecX& z, 
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
        const double tau_squared = tau.dot(tau);
        const double tau_norm = std::sqrt(tau_squared);
        f += tau_norm;

        if (tau_norm > Torque::SQUARE_ROOT_THRESHOLD) {
            if (compute_derivatives) {
                const MatX& ptau_pz = idPtr_->ptau_pz(i);
                const VecX ptau_square_pz = ptau_pz.transpose() * tau;
                grad_f += ptau_square_pz / tau_norm;

                if (compute_hessian) {
                    hess_f += ptau_pz.transpose() * ptau_pz / tau_norm;
                    for (int j = 0; j < tau.size(); j++) {
                        const MatX& ptau_pz_pz = idPtr_->ptau_pz_pz(j, i);
                        hess_f += ptau_pz_pz * tau(j) / tau_norm;
                    }
                    hess_f -= ptau_square_pz * ptau_square_pz.transpose() / std::pow(tau_squared, 1.5);
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