#include "InverseDynamics.h"

namespace RAPTOR {

InverseDynamics::InverseDynamics(const Model& model_input, 
                                 const std::shared_ptr<Trajectories>& trajPtr_input) :
    trajPtr_(trajPtr_input) {   
    N = trajPtr_->N;

    modelPtr_ = std::make_shared<Model>(model_input);
    dataPtr_ = std::make_shared<Data>(model_input);

    tau.resize(1, N);
    ptau_pz.resize(1, N);
    ptau_pz_pz.resize(modelPtr_->nv, N);

    for (int i = 0; i < N; i++) {
        tau(i) = VecX::Zero(modelPtr_->nv);
        ptau_pz(i) = MatX::Zero(modelPtr_->nv, trajPtr_->varLength);

        for (int j = 0; j < modelPtr_->nv; j++) {
            ptau_pz_pz(j, i) = MatX::Zero(trajPtr_->varLength, trajPtr_->varLength);
        }
    }

    prnea_pq = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    prnea_pv = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    prnea_pa = MatX::Zero(modelPtr_->nv, modelPtr_->nv);

    ptau2_pq = Ten3(modelPtr_->nv, modelPtr_->nv, modelPtr_->nv);
    ptau2_pv = Ten3(modelPtr_->nv, modelPtr_->nv, modelPtr_->nv);
    ptau2_pqpv = Ten3(modelPtr_->nv, modelPtr_->nv, modelPtr_->nv);
    ptau2_papq = Ten3(modelPtr_->nv, modelPtr_->nv, modelPtr_->nv);
    ptau2_pq.setZero();
    ptau2_pv.setZero();
    ptau2_pqpv.setZero();
    ptau2_papq.setZero();
}

bool InverseDynamics::is_computed(const VecX& z, bool compute_derivatives, bool compute_hessian) {
    if (compute_hessian && !compute_derivatives) {
        throw std::invalid_argument("compute_derivatives needs to be true when compute_hessian is true.");
        return false;
    }
                            
    if (!Utils::ifTwoVectorEqual(current_z, z, 0)) {
        current_z = z;
        if_compute_derivatives = compute_derivatives;
        if_compute_hessian = compute_hessian;
        return false;
    }

    if (compute_derivatives != if_compute_derivatives) {
        current_z = z;
        if_compute_derivatives = compute_derivatives;
        if_compute_hessian = compute_hessian;
        return false;
    }

    if (compute_hessian != if_compute_hessian) {
        current_z = z;
        if_compute_derivatives = compute_derivatives;
        if_compute_hessian = compute_hessian;
        return false;
    }

    // current_z = z;  
    if_compute_derivatives = compute_derivatives;
    if_compute_hessian = compute_hessian;
    return true;
}

void InverseDynamics::compute(const VecX& z,
                              bool compute_derivatives,
                              bool compute_hessian) {
    if (trajPtr_ == nullptr) {
        throw std::runtime_error("trajPtr_ is not defined yet!");
    }   

    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);                    

    for (int i = 0; i < N; i++) {
        const VecX& q = trajPtr_->q(i);
        const VecX& v = trajPtr_->q_d(i);
        const VecX& a = trajPtr_->q_dd(i);

        const MatX& pq_pz = trajPtr_->pq_pz(i);
        const MatX& pv_pz = trajPtr_->pq_d_pz(i);
        const MatX& pa_pz = trajPtr_->pq_dd_pz(i);

        const Eigen::Array<MatX, Eigen::Dynamic, 1>& pq_pz_pz = trajPtr_->pq_pz_pz.col(i);
        const Eigen::Array<MatX, Eigen::Dynamic, 1>& pv_pz_pz = trajPtr_->pq_d_pz_pz.col(i);
        const Eigen::Array<MatX, Eigen::Dynamic, 1>& pa_pz_pz = trajPtr_->pq_dd_pz_pz.col(i);

        if (compute_hessian) {
            pinocchio::ComputeRNEASecondOrderDerivatives(*modelPtr_, *dataPtr_, 
                                                         q, v, a,
                                                         ptau2_pq, ptau2_pv, ptau2_pqpv, ptau2_papq);
        }
        
        if (compute_derivatives || compute_hessian) {
            pinocchio::computeRNEADerivatives(*modelPtr_, *dataPtr_, 
                                              q, v, a, 
                                              prnea_pq, prnea_pv, prnea_pa);
        }
        else {
            pinocchio::rnea(*modelPtr_, *dataPtr_, q, v, a);
        }
        
        // adjust with damping force
        tau(i) = dataPtr_->tau + 
                 modelPtr_->damping.cwiseProduct(v);
        
        // adjust with static friction
        for (int j = 0; j < tau(i).size(); j++) {
            tau(i)(j) += modelPtr_->friction(j) * Utils::sign(v(j));
        }
                
        if (compute_derivatives) {
            // prnea_pa is just the inertia matrix.
            // pinocchio only computes the upper triangle part of it.
            for (int mi = 0; mi < prnea_pa.rows(); mi++) {
                for (int mj = 0; mj <= mi; mj++) {
                    prnea_pa(mi, mj) = prnea_pa(mj, mi);
                }
            }

            // pinocchio rnea does not take damping into account
            for (int mi = 0; mi < prnea_pa.rows(); mi++) {
                prnea_pv(mi, mi) += modelPtr_->damping(mi);
            }

            ptau_pz(i) = prnea_pq * pq_pz + 
                         prnea_pv * pv_pz + 
                         prnea_pa * pa_pz;

            if (compute_hessian) {
                for (int j = 0; j < tau(i).size(); j++) {
                    const MatX ptau2_pq_local = chipFromTensor3x(ptau2_pq, j, 0);
                    const MatX ptau2_pv_local = chipFromTensor3x(ptau2_pv, j, 0);
                    const MatX ptau2_pqpv_local = chipFromTensor3x(ptau2_pqpv, j, 0);
                    const MatX ptau2_papq_local = chipFromTensor3x(ptau2_papq, j, 0);

                    // p( ptau / pq * pq / pz ) / pz
                    MatX term1 = (ptau2_pq_local * pq_pz + 
                                  ptau2_pqpv_local * pv_pz +
                                  ptau2_papq_local.transpose() * pa_pz).transpose() * pq_pz;
                    for (int k = 0; k < modelPtr_->nv; k++) {
                        term1 += prnea_pq(j, k) * pq_pz_pz(k);
                    }

                    // p( ptau / pv ) / pz * pv / pz
                    MatX term2 = (ptau2_pqpv_local.transpose() * pq_pz + 
                                  ptau2_pv_local * pv_pz).transpose() * pv_pz;
                    for (int k = 0; k < modelPtr_->nv; k++) {
                        term2 += prnea_pv(j, k) * pv_pz_pz(k);
                    }

                    // p( ptau / pa ) / pz * pa / pz
                    MatX term3 = (ptau2_papq_local * pq_pz).transpose() * pa_pz;
                    for (int k = 0; k < modelPtr_->nv; k++) {
                        term3 += prnea_pa(j, k) * pa_pz_pz(k);
                    }

                    // add together
                    ptau_pz_pz(j, i) = term1 + term2 + term3;
                }
            }
        }
    }
}

Eigen::MatrixXd InverseDynamics::chipFromTensor3x(const Ten3& tensor3x, 
                                                  const Eigen::Index offset, 
                                                  const Eigen::Index dim) {
    const auto& dims = tensor3x.dimensions();
    Eigen::Index dim0 = 0, dim1 = 0;
    if (dim == 0) {
        dim0 = dims[1];
        dim1 = dims[2];
    }
    else if (dim == 1) {
        dim0 = dims[0];
        dim1 = dims[2];
    }
    else if (dim == 2) {
        dim0 = dims[0];
        dim1 = dims[1];
    }
    else {
        throw std::invalid_argument("dim should be 0, 1, or 2.");
    }

    Eigen::MatrixXd result(dim0, dim1);

    if (dim == 0) {
        for (Eigen::Index i = 0; i < dim0; i++) {
            for (Eigen::Index j = 0; j < dim1; j++) {
                result(i, j) = tensor3x(offset, i, j);
            }
        }
    }
    else if (dim == 1) {
        for (Eigen::Index i = 0; i < dim0; i++) {
            for (Eigen::Index j = 0; j < dim1; j++) {
                result(i, j) = tensor3x(i, offset, j);
            }
        }
    }
    else if (dim == 2) {
        for (Eigen::Index i = 0; i < dim0; i++) {
            for (Eigen::Index j = 0; j < dim1; j++) {
                result(i, j) = tensor3x(i, j, offset);
            }
        }
    }

    return result;
}

}; // namespace RAPTOR

