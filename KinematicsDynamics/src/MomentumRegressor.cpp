#include "MomentumRegressor.h"

namespace RAPTOR {

void MomentumRegressor::compute(const VecX& z,
                                bool compute_derivatives,
                                bool compute_hessian) {
    if (trajPtr_ == nullptr) {
        throw std::runtime_error("trajPtr_ is not defined yet!");
    }   

    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    if (compute_hessian) {
        throw std::invalid_argument("MomentumRegressor: Hessian not implemented yet");
    }

    int i = 0;
    #pragma omp parallel for shared(trajPtr_, modelPtr_, jtype, Xtree, Y, pY_pz, tau, ptau_pz) private(i) schedule(dynamic, 1)
    for (i = 0; i < N; i++) {
        const VecX& q = trajPtr_->q(i);
        const VecX& q_d = trajPtr_->q_d(i);

        const MatX& pq_pz = trajPtr_->pq_pz(i);
        const MatX& pq_d_pz = trajPtr_->pq_d_pz(i);

        // below is the extended regressor algorithm from ROAM-Lab
        // refer to https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/dynamics/Hqd_and_CTqd.m

        Vec6 vJ;
        Mat6 XJ, dXJdq;
        MatRegressor bodyRegressor;
        Vec6 Sd;
        MatX pSd_pz(6, trajPtr_->varLength);

        Eigen::Array<Mat6, 1, Eigen::Dynamic> Xup(modelPtr_->nv);
        Eigen::Array<Mat6, 1, Eigen::Dynamic> dXupdq(modelPtr_->nv);
        Eigen::Array<Vec6, 1, Eigen::Dynamic> S(modelPtr_->nv);
        Eigen::Array<Vec6, 1, Eigen::Dynamic> v(modelPtr_->nv);
        Eigen::Array<MatX, 1, Eigen::Dynamic> pv_pz(modelPtr_->nv);
        Eigen::Array<MatRegressor, 1, Eigen::Dynamic> pbodyRegressor_pz(trajPtr_->varLength);

        // forward pass
        for (int j = 0; j < modelPtr_->nv; j++) {
            const int pinocchio_joint_id = j + 1; // the first joint in pinocchio is the root joint
            const int parent_id = modelPtr_->parents[pinocchio_joint_id] - 1;

            if (compute_derivatives) {
                Spatial::jcalc(XJ, dXJdq, S(j), jtype(j), q(j));
            }
            else {
                Spatial::jcalc(XJ, S(j), jtype(j), q(j));
            }

            Xup(j) = XJ * Xtree(j);

            if (compute_derivatives) {
                dXupdq(j) = dXJdq * Xtree(j);
            }

            if (parent_id > -1) {
                v(j) = Xup(j) * v(parent_id) + S(j) * q_d(j);

                if (compute_derivatives) {
                    // compute pv_pz
                    pv_pz(j) = Xup(j) * pv_pz(parent_id);

                    if (j < trajPtr_->Nact) {
                        // deal with S(j) * q_d(j);
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k) != 0) {
                                pv_pz(j).row(k) += S(j)(k) * pq_d_pz.row(j);
                            }
                        }

                        // deal with dXupdq(j) * v(parent_id)
                        Vec6 dXupdq_a_parent_id = dXupdq(j) * v(parent_id);
                        for (int k = 0; k < 6; k++) {
                            pv_pz(j).row(k) += dXupdq_a_parent_id(k) * pq_pz.row(j);
                        }
                    }
                }
            }
            else {
                v(j) = S(j) * q_d(j);

                if (compute_derivatives) {// compute pv_pz
                    pv_pz(j) = MatX::Zero(6, trajPtr_->varLength);

                    if (j < trajPtr_->Nact) {
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k) != 0) {
                                pv_pz(j).row(k) = S(j)(k) * pq_d_pz.row(j);
                            }
                        }
                    }
                }
            }

            // compute local bodyRegressor
            const double& v1 = v(j)(0);
            const double& v2 = v(j)(1);
            const double& v3 = v(j)(2);
            const double& v4 = v(j)(3);
            const double& v5 = v(j)(4);
            const double& v6 = v(j)(5);

            bodyRegressor << 
                 0,     0,   v6,   -v5,    v1,    v2,     0,   v3,    0,    0,
                 0,   -v6,    0,    v4,     0,    v1,    v2,    0,   v3,    0,
                 0,    v5,  -v4,     0,     0,     0,     0,   v1,   v2,   v3,
                v4,     0,  -v3,    v2,     0,     0,     0,    0,    0,    0,
                v5,    v3,    0,   -v1,     0,     0,     0,    0,    0,    0,
                v6,   -v2,   v1,     0,     0,     0,     0,    0,    0,    0;
                 
            if (compute_derivatives) {
                for (int k = 0; k < trajPtr_->varLength; k++) {
                    const double& pv1 = pv_pz(j)(0, k);
                    const double& pv2 = pv_pz(j)(1, k);
                    const double& pv3 = pv_pz(j)(2, k);
                    const double& pv4 = pv_pz(j)(3, k);
                    const double& pv5 = pv_pz(j)(4, k);
                    const double& pv6 = pv_pz(j)(5, k);

                    pbodyRegressor_pz(k) <<
                        0,       0,   pv6,   -pv5,   pv1,   pv2,     0,  pv3,    0,    0,
                        0,    -pv6,     0,    pv4,     0,   pv1,   pv2,    0,  pv3,    0,
                        0,     pv5,  -pv4,      0,     0,     0,     0,  pv1,  pv2,  pv3,
                        pv4,     0,  -pv3,    pv2,     0,     0,     0,    0,    0,    0,
                        pv5,   pv3,     0,   -pv1,     0,     0,     0,    0,    0,    0,
                        pv6,  -pv2,   pv1,      0,     0,     0,     0,    0,    0,    0;
                }
            }

            // reverse substitution
            for (int h = j; h >= 0; h = modelPtr_->parents[h + 1] - 1) {
                // regressor for the system momentum
                Y.block(i * modelPtr_->nv + h, j * 10, 1, 10) = 
                    S(h).transpose() * bodyRegressor;

                if (compute_derivatives) {
                    for (int k = 0; k < trajPtr_->varLength; k++) {
                        pY_pz(k).block(i * modelPtr_->nv + h, j * 10, 1, 10) = 
                            S(h).transpose() * pbodyRegressor_pz(k);
                    }
                }

                Sd = Spatial::crm(v(h)) * S(h);

                if (compute_derivatives) {
                    pSd_pz.row(0) = S(h)(2) * pv_pz(h).row(1) - S(h)(1) * pv_pz(h).row(2);
                    pSd_pz.row(1) = S(h)(0) * pv_pz(h).row(2) - S(h)(2) * pv_pz(h).row(0);
                    pSd_pz.row(2) = S(h)(1) * pv_pz(h).row(0) - S(h)(0) * pv_pz(h).row(1);
                    pSd_pz.row(3) = S(h)(2) * pv_pz(h).row(4) - S(h)(1) * pv_pz(h).row(5) - S(h)(4) * pv_pz(h).row(2) + S(h)(5) * pv_pz(h).row(1);
                    pSd_pz.row(4) = S(h)(0) * pv_pz(h).row(5) - S(h)(2) * pv_pz(h).row(3) + S(h)(3) * pv_pz(h).row(2) - S(h)(5) * pv_pz(h).row(0);
                    pSd_pz.row(5) = S(h)(1) * pv_pz(h).row(3) - S(h)(0) * pv_pz(h).row(4) - S(h)(3) * pv_pz(h).row(1) + S(h)(4) * pv_pz(h).row(0);
                }

                // regressor for C^T(q, v) * v
                Y_CTv.block(i * modelPtr_->nv + h, j * 10, 1, 10) = 
                    Sd.transpose() * bodyRegressor;

                if (compute_derivatives) {
                    for (int k = 0; k < trajPtr_->varLength; k++) {
                        pY_CTv_pz(k).block(i * modelPtr_->nv + h, j * 10, 1, 10) = 
                            pSd_pz.col(k).transpose() * bodyRegressor +
                            Sd.transpose() * pbodyRegressor_pz(k);
                    }
                }

                if (modelPtr_->parents[h + 1] > 0) {
                    if (compute_derivatives) {
                        MatRegressor dXupdq_Yfull = dXupdq(h).transpose() * bodyRegressor;
                        for (int k = 0; k < trajPtr_->varLength; k++) {
                            pbodyRegressor_pz(k) = 
                                Xup(h).transpose() * pbodyRegressor_pz(k) + dXupdq_Yfull * pq_pz(h, k);
                        }
                    }

                    bodyRegressor = Xup(h).transpose() * bodyRegressor;
                }
            }
        }

        // compute system momentum
        tau(i) = Y.middleRows(i * modelPtr_->nv, modelPtr_->nv) * phi;

        if (compute_derivatives) {
            for (int k = 0; k < trajPtr_->varLength; k++) {
                ptau_pz(i).col(k) = pY_pz(k).middleRows(i * modelPtr_->nv, modelPtr_->nv) * phi;
            }
        }
    }
}

}; // namespace RAPTOR