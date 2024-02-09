#include "CustomizedInverseDynamics.h"

namespace IDTO {

CustomizedInverseDynamics::CustomizedInverseDynamics(const Model& model_input, 
                                                     const Eigen::VectorXi& jtype_input,
                                                     const std::shared_ptr<Trajectories>& trajPtr_input) {   
    jtype = jtype_input;
    trajPtr_ = trajPtr_input;
    N = trajPtr_->N;

    modelPtr_ = std::make_shared<Model>(model_input);
    dataPtr_ = std::make_shared<Data>(model_input);

    Xtree.resize(1, modelPtr_->nv);
    I.resize(1, modelPtr_->nv);

    for (int i = 0; i < modelPtr_->nv; i++) {
        const int pinocchio_joint_id = i + 1; // the first joint in pinocchio is the root joint

        // plux in Roy Featherstone's code (transformation matrix from parent to child)
        Xtree(i) = plux(modelPtr_->jointPlacements[pinocchio_joint_id].rotation().transpose(), 
                        modelPtr_->jointPlacements[pinocchio_joint_id].translation());
        
        // mcI in Roy Featherstone's code (parallel axis theorem)
        const MatX CC = skew(modelPtr_->inertias[pinocchio_joint_id].lever());
        const double mm = modelPtr_->inertias[pinocchio_joint_id].mass();
        const MatX II = modelPtr_->inertias[pinocchio_joint_id].inertia().matrix();
        I(i) << mm * CC * CC.transpose() + II, mm * CC,
                mm * CC.transpose(),           mm * MatX::Identity(3, 3);
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
    f.resize(1, modelPtr_->nv);
    pf_pz.resize(1, modelPtr_->nv);

    tau.resize(1, N);
    ptau_pz.resize(1, N);
    
    contactForce.resize(1, N);  
    pcontactForce_pz.resize(1, N);
}

void CustomizedInverseDynamics::compute(const VecX& z,
                                        bool compute_derivatives) {
    if (trajPtr_ == nullptr) {
        throw std::runtime_error("trajPtr_ is not defined yet!");
    }   

    if (is_computed(z, compute_derivatives)) {
        return;
    }

    trajPtr_->compute(z, compute_derivatives);                            

    for (int i = 0; i < N - 6; i++) {
        const VecX& q = trajPtr_->q(i);
        const VecX& q_d = trajPtr_->q_d(i);
        const VecX& q_dd = trajPtr_->q_dd(i);

        const MatX& pq_pz = trajPtr_->pq_pz(i);
        const MatX& pq_d_pz = trajPtr_->pq_d_pz(i);
        const MatX& pq_dd_pz = trajPtr_->pq_dd_pz(i);

        tau(i) = VecX::Zero(trajPtr_->Nact);

        // below is the original Roy Featherstone's inverse dynamics algorithm
        // refer to https://royfeatherstone.org/spatial/v2/index.html#ID

        // forward pass
        Vec6 vJ;
        Mat6 XJ, dXJdq;
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

                    if (j < trajPtr_->Nact) {
                        // deal with crm_v_j * vJ;
                        pa_pz(j).row(0) += vJ(2) * pv_pz(j).row(1) - vJ(1) * pv_pz(j).row(2);
                        pa_pz(j).row(1) += vJ(0) * pv_pz(j).row(2) - vJ(2) * pv_pz(j).row(0);
                        pa_pz(j).row(2) += vJ(1) * pv_pz(j).row(0) - vJ(0) * pv_pz(j).row(1);
                        pa_pz(j).row(3) += vJ(2) * pv_pz(j).row(4) - vJ(1) * pv_pz(j).row(5) - vJ(4) * pv_pz(j).row(2) + vJ(5) * pv_pz(j).row(1);
                        pa_pz(j).row(4) += vJ(0) * pv_pz(j).row(5) - vJ(2) * pv_pz(j).row(3) + vJ(3) * pv_pz(j).row(2) - vJ(5) * pv_pz(j).row(0);
                        pa_pz(j).row(5) += vJ(1) * pv_pz(j).row(3) - vJ(0) * pv_pz(j).row(4) - vJ(3) * pv_pz(j).row(1) + vJ(4) * pv_pz(j).row(0);
                        
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
                    pv_pz(j) = MatX::Zero(6, z.size());

                    if (j < trajPtr_->Nact) {
                        for (int k = 0; k < S(j).size(); k++) {
                            if (S(j)(k) != 0) {
                                pv_pz(j).row(k) = S(j)(k) * pq_d_pz.row(j);
                            }
                        }
                    }

                    // compute pa_pz
                    pa_pz(j) = MatX::Zero(6, z.size());

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

            f(j) = I(j) * a(j) + crf(v(j)) * I(j) * v(j);
        }

        // backward pass
        for (int j = modelPtr_->nv - 1; j >= 0; j--) {
            const int pinocchio_joint_id = j + 1; // the first joint in pinocchio is the root joint
            const int parent_id = modelPtr_->parents[pinocchio_joint_id] - 1;

            if (j < trajPtr_->Nact) {
                tau(i)(j) = S(j).transpose() * f(j) + 
                            modelPtr_->rotorInertia(j) * q_dd(j) +
                            modelPtr_->damping(j) * q_d(j) +
                            modelPtr_->friction(j) * sign(q_d(j));
            }
            else {
                contactForce(i) = f(j); // record the contact force at the last joint
            }

            if (parent_id > -1) {
                f(parent_id) += Xup(j).transpose() * f(j);
            }
        }
    }
}

}; // namespace IDTO

