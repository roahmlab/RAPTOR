#include "DigitModifiedCustomizedConstraints.h"

namespace IDTO {
namespace DigitModified {

DigitModifiedCustomizedConstraints::DigitModifiedCustomizedConstraints(const Model& model_input,
                                                                       const Eigen::VectorXi& jtype_input,
                                                                       std::shared_ptr<Trajectories>& trajPtr_input,
                                                                       std::shared_ptr<DynamicsConstraints>& dcPtr_input,
                                                                       const GaitParameters& gp_input) : 
    jtype(jtype_input),
    trajPtr_(trajPtr_input),
    dcPtr_(dcPtr_input),
    gp(gp_input) {
    modelPtr_ = std::make_unique<Model>(model_input);
    fkhofPtr_ = std::make_unique<ForwardKinematicsHighOrderDerivative>();

    if (modelPtr_->existJointName("right_toe_roll")) {
        swingfoot_id = modelPtr_->getJointId("right_toe_roll");
    }
    else {
        throw std::runtime_error("Can not find joint: right_toe_roll");
    }

    // This is right foot end transform
    // This only applies when stance foot is left foot!!!
    swingfoot_endT.R << 0,             -1, 0,
                        0.5,           0,  -sin(M_PI / 3),
                        sin(M_PI / 3), 0,  0.5;
    swingfoot_endT.p << 0, 0, 0;

    jointTJ = MatX::Zero(6, modelPtr_->nv);
    q = MatX::Zero(modelPtr_->nv, trajPtr_->N);
    swingfoot_xyzrpy = MatX::Zero(6, trajPtr_->N);
    pq_pz.resize(1, trajPtr_->N);
    pswingfoot_xyzrpy_pz.resize(1, trajPtr_->N);

    m = trajPtr_->N * 8 + 4;

    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, trajPtr_->varLength);
}

void DigitModifiedCustomizedConstraints::compute(const VecX& z, 
                                                 bool compute_derivatives,
                                                 bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    if (compute_hessian) {
        throw std::invalid_argument("DigitModifiedCustomizedConstraints does not support hessian computation");
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    // compute full joint trajectories and swing foot forward kinematics
    for (int i = 0; i < trajPtr_->N; i++) {
        VecX qi(modelPtr_->nq);
        dcPtr_->fill_independent_vector(qi, trajPtr_->q(i));
        dcPtr_->setupJointPosition(qi, compute_derivatives);
        q.col(i) = qi;

        fkhofPtr_->fk(jointT, *modelPtr_, jtype, swingfoot_id, 0, qi, swingfoot_endT, startT);

        swingfoot_xyzrpy.col(i) = fkhofPtr_->Transform2xyzrpy(jointT);

        if (compute_derivatives) {
            pq_pz(i).resize(modelPtr_->nv, trajPtr_->varLength);
            // fill in independent joints derivatives directly
            for (int j = 0; j < dcPtr_->numIndependentJoints; j++) {
                int indenpendentJointIndex = dcPtr_->return_independent_joint_index(j);
                pq_pz(i).row(indenpendentJointIndex) = trajPtr_->pq_pz(i).row(j);
            }
            // compute and fill in dependent joints derivatives
            MatX pq_dep_pz = dcPtr_->pq_dep_pq_indep * trajPtr_->pq_pz(i);
            for (int j = 0; j < dcPtr_->numDependentJoints; j++) {
                int denpendentJointIndex = dcPtr_->return_dependent_joint_index(j);
                pq_pz(i).row(denpendentJointIndex) = pq_dep_pz.row(j);
            }

            fkhofPtr_->fk_jacobian(dTdq, *modelPtr_, jtype, swingfoot_id, 0, qi, swingfoot_endT, startT);
            fkhofPtr_->Transform2xyzrpyJacobian(jointTJ, jointT, dTdq);
            pswingfoot_xyzrpy_pz(i) = jointTJ * pq_pz(i);
        }
    }

    // (1) swingfoot height always larger than 0
    //               height equal to 0 at the beginning and at the end
    //               height higher than the desired value in the middle
    VecX g1 = swingfoot_xyzrpy.row(2).transpose();

    // (2) swingfoot always flat and points forward
    VecX g2 = swingfoot_xyzrpy.row(3).transpose(); // swingfoot roll
    VecX g3 = swingfoot_xyzrpy.row(4).transpose(); // swingfoot pitch
    VecX g4 = swingfoot_xyzrpy.row(5).transpose().array() + M_PI / 2; // swingfoot yaw

    // (3) swingfoot xy equal to desired value at the beginning and at the end
    VecX g5(4);
    g5 << swingfoot_xyzrpy.topLeftCorner(2, 1), swingfoot_xyzrpy.topRightCorner(2, 1);

    // (4) torso height always larger than 1 meter
    //           roll and pitch always close to 0
    //           yaw always close to 0 when walking forward
    VecX g6 = q.row(2).transpose(); // torso height
    VecX g7 = q.row(3).transpose(); // torso roll
    VecX g8 = q.row(4).transpose(); // torso pitch
    VecX g9 = q.row(5).transpose().array() + M_PI / 2; // torso yaw

    g << g1, g2, g3, g4, g5, g6, g7, g8, g9;

    if (compute_derivatives) {
        int iter = 0;
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(i).row(2);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(i).row(3);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(i).row(4);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(i).row(5);
        }
        pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(0).row(0);
        pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(0).row(1);
        pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(trajPtr_->N - 1).row(0);
        pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(trajPtr_->N - 1).row(1);
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pq_pz(i).row(2);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pq_pz(i).row(3);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pq_pz(i).row(4);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pq_pz(i).row(5);
        }
    }
}

void DigitModifiedCustomizedConstraints::compute_bounds() {
    // (1) swingfoot height always larger than 0
    //               height equal to 0 at the beginning and at the end
    //               height higher than the desired value in the middle
    VecX g1_lb = VecX::Zero(trajPtr_->N);
    VecX g1_ub = VecX::Constant(trajPtr_->N, 1e19);
    g1_lb(trajPtr_->N / 2) = gp.swingfoot_midstep_z_des;
    g1_ub(0) = 0;
    g1_ub(trajPtr_->N - 1) = 0;

    // (2) swingfoot always flat and points forward
    VecX g2_lb = VecX::Zero(trajPtr_->N);
    VecX g2_ub = VecX::Zero(trajPtr_->N);
    VecX g3_lb = VecX::Zero(trajPtr_->N);
    VecX g3_ub = VecX::Zero(trajPtr_->N);
    VecX g4_lb = VecX::Zero(trajPtr_->N);
    VecX g4_ub = VecX::Zero(trajPtr_->N);

    // (3) swingfoot xy equal to desired value at the beginning and at the end
    VecX g5_lb(4);
    VecX g5_ub(4);
    g5_lb << gp.swingfoot_begin_x_des, gp.swingfoot_begin_y_des, gp.swingfoot_end_x_des, gp.swingfoot_end_y_des;
    g5_ub << gp.swingfoot_begin_x_des, gp.swingfoot_begin_y_des, gp.swingfoot_end_x_des, gp.swingfoot_end_y_des;

    // (4) torso height always larger than 1 meter
    //           roll and pitch always close to 0
    //           yaw always close to 0 when walking forward
    VecX g6_lb = VecX::Constant(trajPtr_->N, 1);
    VecX g6_ub = VecX::Constant(trajPtr_->N, 1e19);
    VecX g7_lb = VecX::Constant(trajPtr_->N, -gp.eps_torso_angle);
    VecX g7_ub = VecX::Constant(trajPtr_->N, gp.eps_torso_angle);
    VecX g8_lb = VecX::Constant(trajPtr_->N, -gp.eps_torso_angle);
    VecX g8_ub = VecX::Constant(trajPtr_->N, gp.eps_torso_angle);
    VecX g9_lb = VecX::Constant(trajPtr_->N, -gp.eps_torso_angle);
    VecX g9_ub = VecX::Constant(trajPtr_->N, gp.eps_torso_angle);

    g_lb << g1_lb, g2_lb, g3_lb, g4_lb, g5_lb, g6_lb, g7_lb, g8_lb, g9_lb;
    g_ub << g1_ub, g2_ub, g3_ub, g4_ub, g5_ub, g6_ub, g7_ub, g8_ub, g9_ub;
}

}; // namespace DigitModified
}; // namespace IDTO