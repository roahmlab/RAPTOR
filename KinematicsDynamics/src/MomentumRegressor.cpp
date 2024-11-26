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
        throw std::invalid_argument("CustomizedInverseDynamics: Hessian not implemented yet");
    }

    int i = 0;
    #pragma omp parallel for shared(trajPtr_, modelPtr_, jtype, Xtree, I, Y, pY_pz, tau, ptau_pz) private(i) schedule(dynamic, 1)
    for (i = 0; i < N; i++) {
        const VecX& q = trajPtr_->q(i);
        const VecX& q_d = trajPtr_->q_d(i);

        const MatX& pq_pz = trajPtr_->pq_pz(i);
        const MatX& pq_d_pz = trajPtr_->pq_d_pz(i);

        if (compute_derivatives) {
            ptau_pz(i) = MatX::Zero(trajPtr_->Nact, trajPtr_->varLength);
        }

        // below is the original Roy Featherstone's inverse dynamics algorithm
        // refer to https://royfeatherstone.org/spatial/v2/index.html#ID

        // forward pass
        Mat6 XJ, dXJdq;
        Eigen::Array<Mat6, 1, Eigen::Dynamic> Xup;
        Eigen::Array<Mat6, 1, Eigen::Dynamic> dXupdq;
        Eigen::Array<Vec6, 1, Eigen::Dynamic> S;
        Eigen::Array<Vec6, 1, Eigen::Dynamic> v;
        Eigen::Array<MatX, 1, Eigen::Dynamic> pa_pz;
        Xup.resize(1, modelPtr_->nv);
        dXupdq.resize(1, modelPtr_->nv);
        S.resize(1, modelPtr_->nv);
        v.resize(1, modelPtr_->nv);
        pa_pz.resize(1, modelPtr_->nv);

        // The following has nothing to do with motor dynamics parameters
        MatX Yfull = MatX::Zero(6 * modelPtr_->nv, numInertialParams);
        Eigen::Array<MatX, 1, Eigen::Dynamic> pYfull_pz;
        pYfull_pz.resize(trajPtr_->varLength);
        for (int i = 0; i < trajPtr_->varLength; i++) {
            pYfull_pz(i) = MatX::Zero(6 * modelPtr_->nv, numInertialParams);
        }

        MatX Ycurrent = MatX::Zero(modelPtr_->nv, numInertialParams);

        if (compute_derivatives) {
            for (int k = 0; k < trajPtr_->varLength; k++) {
                pYfull_pz(k).setZero();
            }
        }

        for (int j = 0; j < modelPtr_->nv; j++) {
            const int pinocchio_joint_id = j + 1; // the first joint in pinocchio is the root joint
            const int parent_id = modelPtr_->parents[pinocchio_joint_id] - 1;

            if (compute_derivatives) {
                jcalc(XJ, dXJdq, S(j), jtype(j), q(j));
            }
            else {
                jcalc(XJ, S(j), jtype(j), q(j));
            }

            Xup(j) = XJ * Xtree(j);

            if (compute_derivatives) {
                dXupdq(j) = dXJdq * Xtree(j);
            }

            if (parent_id > -1) {
                v(j) = Xup(j) * v(parent_id) + S(j) * q_d(j);

                if (compute_derivatives) {
                    // compute pa_pz
                    pa_pz(j) = Xup(j) * pa_pz(parent_id);

                    if (j < trajPtr_->Nact) {
                        // deal with S(j) * q_d(j);
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k) != 0) {
                                pa_pz(j).row(k) += S(j)(k) * pq_d_pz.row(j);
                            }
                        }

                        // deal with dXupdq(j) * v(parent_id)
                        Vec6 dXupdq_a_parent_id = dXupdq(j) * v(parent_id);
                        for (int k = 0; k < 6; k++) {
                            pa_pz(j).row(k) += dXupdq_a_parent_id(k) * pq_pz.row(j);
                        }
                    }
                }
            }
            else {
                v(j) = S(j) * q_d(j);

                if (compute_derivatives) {// compute pa_pz
                    pa_pz(j) = MatX::Zero(6, trajPtr_->varLength);

                    if (j < trajPtr_->Nact) {
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k) != 0) {
                                pa_pz(j).row(k) = S(j)(k) * pq_d_pz.row(j);
                            }
                        }
                    }
                }
            }

            const double& v1 = v(j)(0);
            const double& v2 = v(j)(1);
            const double& v3 = v(j)(2);
            const double& v4 = v(j)(3);
            const double& v5 = v(j)(4);
            const double& v6 = v(j)(5);

            Yfull.block(6 * j, 10 * j, 6, 10) << 
                 0,     0,   v6,   -v5,    v1,    v2,     0,   v3,    0,    0,
                 0,   -v6,    0,    v4,     0,    v1,    v2,    0,   v3,    0,
                 0,    v5,  -v4,     0,     0,     0,     0,   v1,   v2,   v3,
                v4,     0,  -v3,    v2,     0,     0,     0,    0,    0,    0,
                v5,    v3,    0,   -v1,     0,     0,     0,    0,    0,    0,
                v6,   -v2,   v1,     0,     0,     0,     0,    0,    0,    0;
     
            if (compute_derivatives) {
                for (int k = 0; k < trajPtr_->varLength; k++) {
                    const double& pv1 = pa_pz(j)(0, k);
                    const double& pv2 = pa_pz(j)(1, k);
                    const double& pv3 = pa_pz(j)(2, k);
                    const double& pv4 = pa_pz(j)(3, k);
                    const double& pv5 = pa_pz(j)(4, k);
                    const double& pv6 = pa_pz(j)(5, k);

                    pYfull_pz(k).block(6 * j, 10 * j, 6, 10) <<
                        0,       0,   pv6,   -pv5,   pv1,   pv2,     0,  pv3,    0,    0,
                        0,    -pv6,     0,    pv4,     0,   pv1,   pv2,    0,  pv3,    0,
                        0,     pv5,  -pv4,      0,     0,     0,     0,  pv1,  pv2,  pv3,
                        pv4,     0,  -pv3,    pv2,     0,     0,     0,    0,    0,    0,
                        pv5,   pv3,     0,   -pv1,     0,     0,     0,    0,    0,    0,
                        pv6,  -pv2,   pv1,      0,     0,     0,     0,    0,    0,    0;
                }
            }
        }

        // backward pass
        for (int j = modelPtr_->nv - 1; j >= 0; j--) {
            const int pinocchio_joint_id = j + 1; // the first joint in pinocchio is the root joint
            const int parent_id = modelPtr_->parents[pinocchio_joint_id] - 1;

            if (parent_id > -1) {
                if (compute_derivatives) {
                    MatX dXupdq_Yfull = dXupdq(j).transpose() * Yfull.middleRows(6 * j, 6);
                    for (int k = 0; k < trajPtr_->varLength; k++) {
                        MatX temp1 = Xup(j).transpose() * pYfull_pz(k).middleRows(6 * j, 6);
                        pYfull_pz(k).middleRows(6 * parent_id, 6) += 
                            temp1 + dXupdq_Yfull * pq_pz(j, k);
                    }
                }

                MatX temp2 = Xup(j).transpose() * Yfull.middleRows(6 * j, 6);
                Yfull.middleRows(6 * parent_id, 6) += temp2;
            }

            Ycurrent.row(j) = S(j).transpose() * Yfull.middleRows(6 * j, 6);

            if (compute_derivatives) {
                for (int k = 0; k < trajPtr_->varLength; k++) {
                    pY_pz(k).block(i * modelPtr_->nv + j, 0, 1, numInertialParams) = S(j).transpose() * pYfull_pz(k).middleRows(6 * j, 6);
                }
            }
        }

        Y.block(i * modelPtr_->nv, 0, modelPtr_->nv, numInertialParams) = Ycurrent;

        // compute system momentum (what is stored in tau is the system momentum) 
        tau(i) = Y.middleRows(i * modelPtr_->nv, modelPtr_->nv) * phi;

        if (compute_derivatives) {
            for (int k = 0; k < trajPtr_->varLength; k++) {
                ptau_pz(i).col(k) = pY_pz(k).middleRows(i * modelPtr_->nv, modelPtr_->nv) * phi;
            }
        }
    }
}

}; // namespace RAPTOR