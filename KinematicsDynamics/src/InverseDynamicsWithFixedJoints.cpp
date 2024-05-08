#include "InverseDynamicsWithFixedJoints.h"

namespace IDTO {

InverseDynamicsWithFixedJoints::InverseDynamicsWithFixedJoints(const Model& model_input, 
                                                               std::shared_ptr<Trajectories>& trajPtr_input) :
    InverseDynamics(model_input, trajPtr_input) {
    lambda.resize(1, N);
    plambda_pz.resize(1, N);
}

void InverseDynamicsWithFixedJoints::compute(const VecX& z,
                                             bool compute_derivatives,
                                             bool compute_hessian) {
    if (trajPtr_ == nullptr) {
        throw std::runtime_error("trajPtr_ is not defined yet!");
    }   

    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    if (compute_hessian) {
        throw std::invalid_argument("InverseDynamicsWithFixedJoints does not support hessian computation");
    }

    if (compute_derivatives) {
        // pinocchio does not store contact wrench gradient information
        // have to compute plambda_pz using numerical gradients
        for (int i = 0; i < N; i++) {
            plambda_pz(i).resize(6, z.size());
        }

        Eigen::Array<Force, 1, Eigen::Dynamic> lambda_plus(1, N);
        Eigen::Array<Force, 1, Eigen::Dynamic> lambda_minus(1, N);
        for (int j = 0; j < z.size(); j++) {
            VecX z_plus = z;
            VecX z_minus = z;
            z_plus(j) += FINITE_DIFFERENCE_DELTA;
            z_minus(j) -= FINITE_DIFFERENCE_DELTA;
            
            compute(z_plus, false);
            lambda_plus = lambda;

            compute(z_minus, false);
            lambda_minus = lambda;

            for (int i = 0; i < N; i++) {
                plambda_pz(i).col(j).head(3) = (lambda_plus(i).linear() - lambda_minus(i).linear()) / 
                                               (2 * FINITE_DIFFERENCE_DELTA);
                plambda_pz(i).col(j).tail(3) = (lambda_plus(i).angular() - lambda_minus(i).angular()) / 
                                               (2 * FINITE_DIFFERENCE_DELTA);
            }
        }
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);                            

    for (int i = 0; i < N; i++) {
        const VecX& q = trajPtr_->q(i);
        const VecX& v = trajPtr_->q_d(i);
        const VecX& a = trajPtr_->q_dd(i);

        if (compute_derivatives) {
            pinocchio::computeRNEADerivatives(*modelPtr_, *dataPtr_, q, v, a, 
                                              prnea_pq, prnea_pv, prnea_pa);
        }

        pinocchio::rnea(*modelPtr_, *dataPtr_, q, v, a);
        
        // adjust with damping force and rotor inertia force
        VecX tau_full = dataPtr_->tau + 
                        modelPtr_->damping.cwiseProduct(v) + 
                        modelPtr_->rotorInertia.cwiseProduct(a);
        tau(i) = tau_full.head(trajPtr_->Nact);

        lambda(i) = dataPtr_->f.back();
                
        if (compute_derivatives) {
            // prnea_pa is just the inertia matrix.
            // pinocchio only computes the upper triangle part of it.
            for (int mi = 0; mi < prnea_pa.rows(); mi++) {
                for (int mj = 0; mj <= mi; mj++) {
                    if (mi == mj) {
                        prnea_pa(mi, mj) += modelPtr_->rotorInertia(mi);
                    }
                    else {
                        prnea_pa(mi, mj) = prnea_pa(mj, mi);
                    }
                }
            }

            // pinocchio rnea does not take damping into account
            for (int mi = 0; mi < prnea_pa.rows(); mi++) {
                prnea_pv(mi, mi) += modelPtr_->damping(mi);
            }

            ptau_pz(i) = prnea_pq.topLeftCorner(trajPtr_->Nact, trajPtr_->Nact) * trajPtr_->pq_pz(i) + 
                         prnea_pv.topLeftCorner(trajPtr_->Nact, trajPtr_->Nact) * trajPtr_->pq_d_pz(i) + 
                         prnea_pa.topLeftCorner(trajPtr_->Nact, trajPtr_->Nact) * trajPtr_->pq_dd_pz(i);
        }
    }
}

}; // namespace IDTO

