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
        const VecX& q_dd = trajPtr_->q_d(i);

        const MatX& pq_pz = trajPtr_->pq_pz(i);
        const MatX& pq_dd_pz = trajPtr_->pq_d_pz(i);

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
        Eigen::Array<Vec6, 1, Eigen::Dynamic> a;
        Eigen::Array<MatX, 1, Eigen::Dynamic> pa_pz;
        Xup.resize(1, modelPtr_->nv);
        dXupdq.resize(1, modelPtr_->nv);
        S.resize(1, modelPtr_->nv);
        a.resize(1, modelPtr_->nv);
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
                a(j) = Xup(j) * a(parent_id) + S(j) * q_dd(j);

                if (compute_derivatives) {
                    // compute pa_pz
                    pa_pz(j) = Xup(j) * pa_pz(parent_id);

                    if (j < trajPtr_->Nact) {
                        // deal with S(j) * q_dd(j);
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k) != 0) {
                                pa_pz(j).row(k) += S(j)(k) * pq_dd_pz.row(j);
                            }
                        }

                        // deal with dXupdq(j) * a(parent_id)
                        Vec6 dXupdq_a_parent_id = dXupdq(j) * a(parent_id);
                        for (int k = 0; k < 6; k++) {
                            pa_pz(j).row(k) += dXupdq_a_parent_id(k) * pq_pz.row(j);
                        }
                    }
                }
            }
            else {
                a(j) = S(j) * q_dd(j);

                if (compute_derivatives) {// compute pa_pz
                    pa_pz(j) = MatX::Zero(6, trajPtr_->varLength);

                    if (j < trajPtr_->Nact) {
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k) != 0) {
                                pa_pz(j).row(k) = S(j)(k) * pq_dd_pz.row(j);
                            }
                        }
                    }
                }
            }

            const double& a1 = a(j)(0);
            const double& a2 = a(j)(1);
            const double& a3 = a(j)(2);
            const double& a4 = a(j)(3);
            const double& a5 = a(j)(4);
            const double& a6 = a(j)(5);

            Yfull.block(6 * j, 10 * j, 6, 10) << 
                 0,     0,   a6,   -a5,    a1,    a2,     0,   a3,    0,    0,
                 0,   -a6,    0,    a4,     0,    a1,    a2,    0,   a3,    0,
                 0,    a5,  -a4,     0,     0,     0,     0,   a1,   a2,   a3,
                a4,     0,  -a3,    a2,     0,     0,     0,    0,    0,    0,
                a5,    a3,    0,   -a1,     0,     0,     0,    0,    0,    0,
                a6,   -a2,   a1,     0,     0,     0,     0,    0,    0,    0;
     
            if (compute_derivatives) {
                for (int k = 0; k < trajPtr_->varLength; k++) {
                    const double& pa1 = pa_pz(j)(0, k);
                    const double& pa2 = pa_pz(j)(1, k);
                    const double& pa3 = pa_pz(j)(2, k);
                    const double& pa4 = pa_pz(j)(3, k);
                    const double& pa5 = pa_pz(j)(4, k);
                    const double& pa6 = pa_pz(j)(5, k);

                    pYfull_pz(k).block(6 * j, 10 * j, 6, 10) <<
                        0,       0,   pa6,   -pa5,   pa1,   pa2,     0,  pa3,    0,    0,
                        0,    -pa6,     0,    pa4,     0,   pa1,   pa2,    0,  pa3,    0,
                        0,     pa5,  -pa4,      0,     0,     0,     0,  pa1,  pa2,  pa3,
                        pa4,     0,  -pa3,    pa2,     0,     0,     0,    0,    0,    0,
                        pa5,   pa3,     0,   -pa1,     0,     0,     0,    0,    0,    0,
                        pa6,  -pa2,   pa1,      0,     0,     0,     0,    0,    0,    0;
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

        // // motor dynamics
        // for (int j = 0; j < modelPtr_->nv; j++) {
        //     Y(i * modelPtr_->nv + j, numInertialParams + 3 * j    ) = Utils::sign(q_d(j));
        //     Y(i * modelPtr_->nv + j, numInertialParams + 3 * j + 1) = q_d(j);
        //     Y(i * modelPtr_->nv + j, numInertialParams + 3 * j + 2) = q_dd(j);

        //     if (compute_derivatives) {
        //         for (int k = 0; k < trajPtr_->varLength; k++) {
        //             pY_pz(k)(i * modelPtr_->nv + j, numInertialParams + 3 * j    ) = 0;
        //             pY_pz(k)(i * modelPtr_->nv + j, numInertialParams + 3 * j + 1) = pq_d_pz(j, k);
        //             pY_pz(k)(i * modelPtr_->nv + j, numInertialParams + 3 * j + 2) = pq_dd_pz(j, k);
        //         }
        //     }
        // }

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