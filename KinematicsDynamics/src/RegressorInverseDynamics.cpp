#include "RegressorInverseDynamics.h"

namespace RAPTOR {

RegressorInverseDynamics::RegressorInverseDynamics(const Model& model_input, 
                                                   const std::shared_ptr<Trajectories>& trajPtr_input,
                                                   Eigen::VectorXi jtype_input) :
    jtype(jtype_input) {
    trajPtr_ = trajPtr_input;
    N = trajPtr_->N;
    modelPtr_ = std::make_shared<Model>(model_input);
    dataPtr_ = std::make_shared<Data>(model_input);

    if (jtype.size() > 0) {
        if (modelPtr_->nv != jtype.size()) {
            std::cerr << "modelPtr_->nv = " << modelPtr_->nv << std::endl;
            std::cerr << "jtype.size() = " << jtype.size() << std::endl;
            throw std::invalid_argument("modelPtr_->nv != jtype.size()");
        }
    }
    else {
        jtype = convertPinocchioJointType(*modelPtr_);
    }

    Xtree.resize(1, modelPtr_->nv);
    I.resize(1, modelPtr_->nv);

    numInertialParams = 10 * modelPtr_->nv;
    numParams = 10 * modelPtr_->nv + // inertial parameters
                3 * modelPtr_->nv;   // motor dynamics parameters
    // numParams = 10 * modelPtr_->nv; // inertial parameters
    phi = VecX::Zero(numParams); 

    for (int i = 0; i < modelPtr_->nv; i++) {
        const int pinocchio_joint_id = i + 1; // the first joint in pinocchio is the root joint

        // plux in Roy Featherstone's code (transformation matrix from parent to child)
        Xtree(i) = Utils::plux(modelPtr_->jointPlacements[pinocchio_joint_id].rotation().transpose(), 
                               modelPtr_->jointPlacements[pinocchio_joint_id].translation());
        
        // function mcI in Roy Featherstone's code (parallel axis theorem)
        const Mat3 CC = Utils::skew(modelPtr_->inertias[pinocchio_joint_id].lever());
        const double mm = modelPtr_->inertias[pinocchio_joint_id].mass();
        const Mat3 II = modelPtr_->inertias[pinocchio_joint_id].inertia().matrix();

        // apply parallel axis theorem
        Mat3 mC = mm * CC;
        I(i) << mm * CC * CC.transpose() + II, mC,
                mC.transpose(),                mm * MatX::Identity(3, 3);

        phi.segment<10>(10 * i) << I(i)(0, 0),      // inertia (in parent joint frame)
                                   I(i)(0, 1),
                                   I(i)(0, 2),
                                   I(i)(1, 1),
                                   I(i)(1, 2),
                                   I(i)(2, 2),
                                   Utils::skew(mC), // center of mass (in parent joint frame)
                                   mm;              // mass

        phi.segment<3>(numInertialParams + 3 * i) << modelPtr_->friction(i),
                                                     modelPtr_->damping(i),
                                                     modelPtr_->rotorInertia(i);
    }

    a_grav << modelPtr_->gravity.angular(),
              modelPtr_->gravity.linear();

    Xup.resize(1, modelPtr_->nv);
    dXupdq.resize(1, modelPtr_->nv);
    S.resize(1, modelPtr_->nv);
    v.resize(1, modelPtr_->nv);
    pv_pz.resize(1, modelPtr_->nv);
    a.resize(1, modelPtr_->nv);
    pa_pz.resize(1, modelPtr_->nv);

    tau.resize(1, N);
    ptau_pz.resize(1, N);

    for (int i = 0; i < N; i++) {
        tau(i) = VecX::Zero(modelPtr_->nv);
        ptau_pz(i) = MatX::Zero(modelPtr_->nv, trajPtr_->varLength);
    }

    Y = MatX::Zero(trajPtr_->N * modelPtr_->nv, numParams);
    pY_pz.resize(trajPtr_->varLength);
    for (int i = 0; i < trajPtr_->varLength; i++) {
        pY_pz(i) = MatX::Zero(trajPtr_->N * modelPtr_->nv, numParams);
    }

    // The following has nothing to do with motor dynamics parameters
    Yfull = MatX::Zero(6 * modelPtr_->nv, numInertialParams);
    pYfull_pz.resize(trajPtr_->varLength);
    for (int i = 0; i < trajPtr_->varLength; i++) {
        pYfull_pz(i) = MatX::Zero(6 * modelPtr_->nv, numInertialParams);
    }

    Ycurrent = MatX::Zero(modelPtr_->nv, numInertialParams);
}

void RegressorInverseDynamics::compute(const VecX& z,
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

    for (int i = 0; i < N; i++) {
        const VecX& q = trajPtr_->q(i);
        const VecX& q_d = trajPtr_->q_d(i);
        const VecX& q_dd = trajPtr_->q_dd(i);

        const MatX& pq_pz = trajPtr_->pq_pz(i);
        const MatX& pq_d_pz = trajPtr_->pq_d_pz(i);
        const MatX& pq_dd_pz = trajPtr_->pq_dd_pz(i);

        if (compute_derivatives) {
            ptau_pz(i) = MatX::Zero(trajPtr_->Nact, trajPtr_->varLength);
        }

        // below is the original Roy Featherstone's inverse dynamics algorithm
        // refer to https://royfeatherstone.org/spatial/v2/index.html#ID

        // forward pass
        Vec6 vJ;
        Mat6 XJ, dXJdq;
        Yfull.setZero();
        Ycurrent.setZero();

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

            vJ = S(j) * q_d(j);
            Xup(j) = XJ * Xtree(j);

            if (compute_derivatives) {
                dXupdq(j) = dXJdq * Xtree(j);
            }

            if (parent_id > -1) {
                v(j) = Xup(j) * v(parent_id) + vJ;
                Mat6 crm_v_j = crm(v(j));
                a(j) = Xup(j) * a(parent_id) + crm_v_j * vJ + S(j) * q_dd(j);

                if (compute_derivatives) {
                    // compute pv_pz
                    pv_pz(j) = Xup(j) * pv_pz(parent_id);

                    if (j < trajPtr_->Nact) {
                        // deal with vJ = S(j) * q_d(j);
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k) != 0) {
                                pv_pz(j).row(k) += S(j)(k) * pq_d_pz.row(j);
                            }
                        }

                        // deal with dXupdq(j) * v(parent_id)
                        Vec6 dXupdq_v_parent_id = dXupdq(j) * v(parent_id);
                        for (int k = 0; k < 6; k++) {
                            pv_pz(j).row(k) += dXupdq_v_parent_id(k) * pq_pz.row(j);
                        }
                    }

                    // compute pa_pz
                    pa_pz(j) = Xup(j) * pa_pz(parent_id);

                    // deal with crm_v_j * vJ;
                    pa_pz(j).row(0) += vJ(2) * pv_pz(j).row(1) - vJ(1) * pv_pz(j).row(2);
                    pa_pz(j).row(1) += vJ(0) * pv_pz(j).row(2) - vJ(2) * pv_pz(j).row(0);
                    pa_pz(j).row(2) += vJ(1) * pv_pz(j).row(0) - vJ(0) * pv_pz(j).row(1);
                    pa_pz(j).row(3) += vJ(2) * pv_pz(j).row(4) - vJ(1) * pv_pz(j).row(5) - vJ(4) * pv_pz(j).row(2) + vJ(5) * pv_pz(j).row(1);
                    pa_pz(j).row(4) += vJ(0) * pv_pz(j).row(5) - vJ(2) * pv_pz(j).row(3) + vJ(3) * pv_pz(j).row(2) - vJ(5) * pv_pz(j).row(0);
                    pa_pz(j).row(5) += vJ(1) * pv_pz(j).row(3) - vJ(0) * pv_pz(j).row(4) - vJ(3) * pv_pz(j).row(1) + vJ(4) * pv_pz(j).row(0);

                    if (j < trajPtr_->Nact) {
                        Vec6 crm_v_j_S_j = crm_v_j * S(j);
                        for (int k = 0; k < 6; k++) {
                            pa_pz(j).row(k) += crm_v_j_S_j(k) * pq_d_pz.row(j);
                        }

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
                v(j) = vJ;
                a(j) = Xup(j) * (-a_grav) + S(j) * q_dd(j);

                if (compute_derivatives) {
                    // compute pv_pz
                    pv_pz(j) = MatX::Zero(6, trajPtr_->varLength);

                    if (j < trajPtr_->Nact) {
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k) != 0) {
                                pv_pz(j).row(k) = S(j)(k) * pq_d_pz.row(j);
                            }
                        }
                    }

                    // compute pa_pz
                    pa_pz(j) = MatX::Zero(6, trajPtr_->varLength);

                    if (j < trajPtr_->Nact) {
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k) != 0) {
                                pa_pz(j).row(k) = S(j)(k) * pq_dd_pz.row(j);
                            }
                        }

                        Vec6 dXupdq_a_grav = dXupdq(j) * (-a_grav);
                        for (int k = 0; k < 6; k++) {
                            pa_pz(j).row(k) += dXupdq_a_grav(k) * pq_pz.row(j);
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

            const double& a1 = a(j)(0);
            const double& a2 = a(j)(1);
            const double& a3 = a(j)(2);
            const double& a4 = a(j)(3);
            const double& a5 = a(j)(4);
            const double& a6 = a(j)(5);

            Yfull.block(6 * j, 10 * j, 6, 10) << 
                a1,    a2 - v1*v3,    a3 + v1*v2, -v2*v3, v2*v2 - v3*v3,  v2*v3,                  0, a6 + v1*v5 - v2*v4, v1*v6 - a5 - v3*v4,                  0,
             v1*v3,    a1 + v2*v3, v3*v3 - v1*v1,     a2,    a3 - v1*v2, -v1*v3, v2*v4 - v1*v5 - a6,                  0, a4 + v2*v6 - v3*v5,                  0,
            -v1*v2, v1*v1 - v2*v2,    a1 - v2*v3,  v1*v2,    a2 + v1*v3,     a3, a5 - v1*v6 + v3*v4, v3*v5 - v2*v6 - a4,                  0,                  0,
                 0,             0,             0,      0,             0,      0,     -v2*v2 - v3*v3,         v1*v2 - a3,         a2 + v1*v3, a4 + v2*v6 - v3*v5,
                 0,             0,             0,      0,             0,      0,         a3 + v1*v2,     -v1*v1 - v3*v3,         v2*v3 - a1, a5 - v1*v6 + v3*v4,
                 0,             0,             0,      0,             0,      0,         v1*v3 - a2,         a1 + v2*v3,     -v1*v1 - v2*v2, a6 + v1*v5 - v2*v4;
                 
            if (compute_derivatives) {
                for (int k = 0; k < trajPtr_->varLength; k++) {
                    const double& pv1 = pv_pz(j)(0, k);
                    const double& pv2 = pv_pz(j)(1, k);
                    const double& pv3 = pv_pz(j)(2, k);
                    const double& pv4 = pv_pz(j)(3, k);
                    const double& pv5 = pv_pz(j)(4, k);
                    const double& pv6 = pv_pz(j)(5, k);

                    const double& pa1 = pa_pz(j)(0, k);
                    const double& pa2 = pa_pz(j)(1, k);
                    const double& pa3 = pa_pz(j)(2, k);
                    const double& pa4 = pa_pz(j)(3, k);
                    const double& pa5 = pa_pz(j)(4, k);
                    const double& pa6 = pa_pz(j)(5, k);

                    pYfull_pz(k).block(6 * j, 10 * j, 6, 10) <<
                                 pa1, pa2 - v1*pv3 - pv1*v3, pa3 + v1*pv2 + pv1*v2, -pv2*v3 - v2*pv3,   2*v2*pv2 - 2*v3*pv3,  pv2*v3 + v2*pv3,                                       0, pa6 + v1*pv5 + pv1*v5 - v2*pv4 - pv2*v4, v1*pv6 + pv1*v6 - pa5 - v3*pv4 - pv3*v4,                                       0,
                     v1*pv3 + pv1*v3, pa1 + v2*pv3 + pv2*v3,   2*v3*pv3 - 2*v1*pv1,              pa2, pa3 - v1*pv2 - pv1*v2, -v1*pv3 - pv1*v3, v2*pv4 + pv2*v4 - v1*pv5 - pv1*v5 - pa6,                                       0, pa4 + v2*pv6 + pv2*v6 - v3*pv5 - pv3*v5,                                       0,
                    -v1*pv2 - pv1*v2,   2*v1*pv1 - 2*v2*pv2, pa1 - v2*pv3 - pv2*v3,  v1*pv2 + pv1*v2, pa2 + v1*pv3 + pv1*v3,              pa3, v3*pv4 + pv3*v4 - v1*pv6 - pv1*v6 + pa5, v3*pv5 + pv3*v5 - v2*pv6 - pv2*v6 - pa4,                                       0,                                       0,
                                   0,                     0,                     0,                0,                     0,                0,                    -2*v2*pv2 - 2*v3*pv3,                   v1*pv2 + pv1*v2 - pa3,                   pa2 + v1*pv3 + pv1*v3, pa4 + v2*pv6 + pv2*v6 - v3*pv5 - pv3*v5,
                                   0,                     0,                     0,                0,                     0,                0,                   pa3 + v1*pv2 + pv1*v2,                    -2*v1*pv1 - 2*v3*pv3,                   v2*pv3 + pv2*v3 - pa1, v3*pv4 + pv3*v4 - v1*pv6 - pv1*v6 + pa5,
                                   0,                     0,                     0,                0,                     0,                0,                   v1*pv3 + pv1*v3 - pa2,                   pa1 + v2*pv3 + pv2*v3,                    -2*v1*pv1 - 2*v2*pv2, pa6 + v1*pv5 + pv1*v5 - v2*pv4 - pv2*v4;
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

        // motor dynamics
        for (int j = 0; j < modelPtr_->nv; j++) {
            Y(i * modelPtr_->nv + j, numInertialParams + 3 * j    ) = Utils::sign(q_d(j));
            Y(i * modelPtr_->nv + j, numInertialParams + 3 * j + 1) = q_d(j);
            Y(i * modelPtr_->nv + j, numInertialParams + 3 * j + 2) = q_dd(j);

            if (compute_derivatives) {
                for (int k = 0; k < trajPtr_->varLength; k++) {
                    pY_pz(k)(i * modelPtr_->nv + j, numInertialParams + 3 * j    ) = 0;
                    pY_pz(k)(i * modelPtr_->nv + j, numInertialParams + 3 * j + 1) = pq_d_pz(j, k);
                    pY_pz(k)(i * modelPtr_->nv + j, numInertialParams + 3 * j + 2) = pq_dd_pz(j, k);
                
                }
            }
        }

        // compute torque
        tau(i) = Y.middleRows(i * modelPtr_->nv, modelPtr_->nv) * phi;

        if (compute_derivatives) {
            for (int k = 0; k < trajPtr_->varLength; k++) {
                ptau_pz(i).col(k) = pY_pz(k).middleRows(i * modelPtr_->nv, modelPtr_->nv) * phi;
            }
        }
    }
}

}; // namespace RAPTOR