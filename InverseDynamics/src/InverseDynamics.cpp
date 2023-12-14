#include "InverseDynamics.h"

namespace IDTO {

InverseDynamics::InverseDynamics(const Model& model_input, 
                                 std::shared_ptr<Trajectories>& trajPtr_input) :
    trajPtr_(trajPtr_input) {   
    N = trajPtr_->N;

    modelPtr_ = std::make_unique<Model>(model_input);
    dataPtr_ = std::make_unique<Data>(model_input);

    tau.resize(1, N);

    ptau_pz.resize(1, N);

    prnea_pq = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    prnea_pv = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    prnea_pa = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
}

void InverseDynamics::compute(const VecX& z,
                              bool compute_derivatives) {
    if (trajPtr_ == nullptr) {
        throw std::runtime_error("trajPtr_ is not defined yet!");
    }   

    if (is_computed(z, compute_derivatives)) {
        return;
    }

    trajPtr_->compute(z, compute_derivatives);                            

    for (int i = 0; i < N; i++) {
        const VecX& q = trajPtr_->q(i);
        const VecX& v = trajPtr_->q_d(i);
        const VecX& a = trajPtr_->q_dd(i);

        if (!compute_derivatives) {
            pinocchio::rnea(*modelPtr_, *dataPtr_, q, v, a);
        }
        else {
            pinocchio::computeRNEADerivatives(*modelPtr_, *dataPtr_, q, v, a, 
                                              prnea_pq, prnea_pv, prnea_pa);
        }
        
        // adjust with damping force and rotor inertia force
        tau(i) = dataPtr_->tau + 
                 modelPtr_->damping.cwiseProduct(v) + 
                 modelPtr_->rotorInertia.cwiseProduct(a);
                
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

            ptau_pz(i) = prnea_pq * trajPtr_->pq_pz(i) + 
                         prnea_pv * trajPtr_->pq_d_pz(i) + 
                         prnea_pa * trajPtr_->pq_dd_pz(i);
        }
    }
}

}; // namespace IDTO

