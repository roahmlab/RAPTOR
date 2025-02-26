#include "CustomizedInverseDynamics.h"

namespace RAPTOR {

CustomizedInverseDynamics::CustomizedInverseDynamics(const Model& model_input, 
                                                     const std::shared_ptr<Trajectories>& trajPtr_input,
                                                     Eigen::VectorXi jtype_input) : 
    jtype(jtype_input) {
    trajPtr_ = trajPtr_input;
    N = trajPtr_->N;

    modelPtr_ = std::make_shared<Model>(model_input);
    dataPtr_ = std::make_shared<Data>(model_input);

    // sanity check on consistency on jtype
    if (jtype.size() > 0) {
        if (modelPtr_->nv != jtype.size()) {
            std::cerr << "model.nv: " << modelPtr_->nv << std::endl;
            std::cerr << "jtype.size(): " << jtype.size() << std::endl;
            throw std::invalid_argument("CustomizedInverseDynamics: total number of joints in jtype not consistent with model!");
        }
    }
    else {
        jtype = convertPinocchioJointType(*modelPtr_);
    }

    active_joints.clear();
    for (int i = 0; i < jtype.size(); i++) {
        if (jtype(i) != 0) {
            active_joints.push_back(i);
        }
    }
    if (active_joints.size() != trajPtr_->Nact) {
        std::cerr << "active_joints: " << active_joints.size() << std::endl;
        std::cerr << "trajPtr_->Nact: " << trajPtr_->Nact << std::endl;
        throw std::invalid_argument("CustomizedInverseDynamics: actuated joints in jtype not consistent with trajectory pointer!");
    }

    // initialize the arrays for kinematics and dynamics
    Xtree.resize(1, modelPtr_->nv);
    I.resize(1, modelPtr_->nv);

    for (int i = 0; i < modelPtr_->nv; i++) {
        const int pinocchio_joint_id = i + 1; // the first joint in pinocchio is the root joint

        // plux in Roy Featherstone's code (transformation matrix from parent to child)
        Xtree(i) = Utils::plux(modelPtr_->jointPlacements[pinocchio_joint_id].rotation().transpose(), 
                               modelPtr_->jointPlacements[pinocchio_joint_id].translation());
        
        // function mcI in Roy Featherstone's code (parallel axis theorem)
        const MatX CC = Utils::skew(modelPtr_->inertias[pinocchio_joint_id].lever());
        const double mm = modelPtr_->inertias[pinocchio_joint_id].mass();
        const MatX II = modelPtr_->inertias[pinocchio_joint_id].inertia().matrix();
        I(i) << mm * CC * CC.transpose() + II, mm * CC,
                mm * CC.transpose(),           mm * MatX::Identity(3, 3);
    }

    a_grav << modelPtr_->gravity.angular(),
              modelPtr_->gravity.linear();

    // allocate memory for the arrays
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
    
    lambda.resize(1, N);  
    plambda_pz.resize(1, N);
}

Eigen::VectorXd CustomizedInverseDynamics::get_full_joints(const VecX& q) const {
    Eigen::VectorXd q_full(modelPtr_->nq);
    q_full.setZero();
    for (int i = 0; i < active_joints.size(); i++) {
        q_full(active_joints[i]) = q(i);
    }
    return q_full;
}

Eigen::MatrixXd CustomizedInverseDynamics::get_full_joints_derivative(const MatX& q) const {
    Eigen::MatrixXd q_full_derivative(modelPtr_->nq, q.cols());
    q_full_derivative.setZero();
    for (int i = 0; i < active_joints.size(); i++) {
        q_full_derivative.row(active_joints[i]) = q.row(i);
    }
    return q_full_derivative;
}

void CustomizedInverseDynamics::compute(const VecX& z,
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
        const VecX q = get_full_joints(trajPtr_->q(i));
        const VecX q_d = get_full_joints(trajPtr_->q_d(i));
        const VecX q_dd = get_full_joints(trajPtr_->q_dd(i));

        const MatX pq_pz = get_full_joints_derivative(trajPtr_->pq_pz(i));
        const MatX pq_d_pz = get_full_joints_derivative(trajPtr_->pq_d_pz(i));
        const MatX pq_dd_pz = get_full_joints_derivative(trajPtr_->pq_dd_pz(i));

        tau(i) = VecX::Zero(trajPtr_->Nact);

        if (compute_derivatives) {
            ptau_pz(i) = MatX::Zero(trajPtr_->Nact, trajPtr_->varLength);
            plambda_pz(i) = MatX::Zero(6, trajPtr_->varLength);
        }

        // below is the original Roy Featherstone's inverse dynamics algorithm
        // refer to https://royfeatherstone.org/spatial/v2/index.html#ID

        // forward pass
        Vec6 vJ;
        Mat6 XJ, dXJdq;
        for (int j = 0; j < modelPtr_->nv; j++) {
            const int pinocchio_joint_id = j + 1; // the first joint in pinocchio is the root joint
            const int parent_id = modelPtr_->parents[pinocchio_joint_id] - 1;

            if (compute_derivatives) {
                Spatial::jcalc(XJ, dXJdq, S(j), jtype(j), q(j));
            }
            else {
                Spatial::jcalc(XJ, S(j), jtype(j), q(j));
            }

            vJ = S(j) * q_d(j);
            Xup(j) = XJ * Xtree(j);

            if (compute_derivatives) {
                dXupdq(j) = dXJdq * Xtree(j);
            }

            if (parent_id > -1) {
                v(j) = Xup(j) * v(parent_id) + vJ;
                Mat6 crm_v_j = Spatial::crm(v(j));
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

            Mat6 crf_v_j = Spatial::crf(v(j));
            Vec6 I_j_v_j = I(j) * v(j);
            f(j) = I(j) * a(j) + crf_v_j * I_j_v_j;

            if (compute_derivatives) {
                // compute pf_pz
                pf_pz(j) = I(j) * pa_pz(j) + crf_v_j * I(j) * pv_pz(j);

                // deal with d(crf_v_j) * I(j) * v(j)
                pf_pz(j).row(0) += I_j_v_j(2) * pv_pz(j).row(1) - I_j_v_j(1) * pv_pz(j).row(2) - I_j_v_j(4) * pv_pz(j).row(5) + I_j_v_j(5) * pv_pz(j).row(4);
                pf_pz(j).row(1) += I_j_v_j(0) * pv_pz(j).row(2) - I_j_v_j(2) * pv_pz(j).row(0) + I_j_v_j(3) * pv_pz(j).row(5) - I_j_v_j(5) * pv_pz(j).row(3);
                pf_pz(j).row(2) += I_j_v_j(1) * pv_pz(j).row(0) - I_j_v_j(0) * pv_pz(j).row(1) - I_j_v_j(3) * pv_pz(j).row(4) + I_j_v_j(4) * pv_pz(j).row(3);
                pf_pz(j).row(3) += I_j_v_j(5) * pv_pz(j).row(1) - I_j_v_j(4) * pv_pz(j).row(2);
                pf_pz(j).row(4) += I_j_v_j(3) * pv_pz(j).row(2) - I_j_v_j(5) * pv_pz(j).row(0);
                pf_pz(j).row(5) += I_j_v_j(4) * pv_pz(j).row(0) - I_j_v_j(3) * pv_pz(j).row(1);
            }
        }

        // backward pass
        for (int j = modelPtr_->nv - 1; j >= 0; j--) {
            const int pinocchio_joint_id = j + 1; // the first joint in pinocchio is the root joint
            const int parent_id = modelPtr_->parents[pinocchio_joint_id] - 1;

            if (j < trajPtr_->Nact) {
                tau(i)(j) = S(j).transpose() * f(j) + 
                            modelPtr_->armature(j) * q_dd(j) +
                            modelPtr_->damping(j) * q_d(j) +
                            modelPtr_->friction(j) * Utils::sign(q_d(j));

                if (compute_derivatives) {
                    ptau_pz(i).row(j) = S(j).transpose() * pf_pz(j) + 
                                        modelPtr_->armature(j) * pq_dd_pz.row(j) +
                                        modelPtr_->damping(j) * pq_d_pz.row(j);
                }
            }
            else {
                lambda(i) = f(j); // record the contact force at the last joint

                if (compute_derivatives) {
                    plambda_pz(i) = pf_pz(j);
                }
            }

            if (parent_id > -1) {
                f(parent_id) += Xup(j).transpose() * f(j);

                if (compute_derivatives) {
                    pf_pz(parent_id) += Xup(j).transpose() * pf_pz(j);

                    // deal with Xup(j).transpose() * f(j)
                    if (j < trajPtr_->Nact) {
                        Vec6 dXupdq_T_f_j = dXupdq(j).transpose() * f(j);
                        for (int k = 0; k < 6; k++) {
                            pf_pz(parent_id).row(k) += dXupdq_T_f_j(k) * pq_pz.row(j);
                        }
                    }
                }
            }
        }
    }
}

}; // namespace RAPTOR

