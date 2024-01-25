#include "DigitModifiedDynamicsConstraints.h"

namespace IDTO {
namespace DigitModified {

DigitModifiedDynamicsConstraints::DigitModifiedDynamicsConstraints(const Model& model_input, 
                                                                   const Eigen::VectorXi& jtype_input,
                                                                   char stanceLeg_input, 
                                                                   const Transform& stance_foot_T_des_input) :
    DynamicsConstraints(model_input, NUM_DEPENDENT_JOINTS),
    jtype(jtype_input),
    stanceLeg(stanceLeg_input) {
    fkhofPtr_ = std::make_unique<ForwardKinematicsHighOrderDerivative>();

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        if (modelPtr_->existJointName(dependentJointNames[i])) {
            dependentJointIds[i] = modelPtr_->getJointId(dependentJointNames[i]) - 1;
        }
        else {
            throw std::runtime_error("Can not find joint: " + dependentJointNames[i]);
        }
    }
    
    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        if (modelPtr_->existJointName(independentJointNames[i])) {
            independentJointIds[i] = modelPtr_->getJointId(independentJointNames[i]) - 1;
        }
        else {
            throw std::runtime_error("Can not find joint: " + dependentJointNames[i]);
        }
    }

    if (stanceLeg == 'L' || stanceLeg == 'l') {
        if (modelPtr_->existJointName("left_toe_roll")) {
            contact_joint_id = modelPtr_->getJointId("left_toe_roll");
        }
        else {
            throw std::runtime_error("Can not find joint: left_toe_roll");
        }

        stance_foot_endT.R << 0,             1, 0,
                              -0.5,          0, sin(M_PI / 3),
                              sin(M_PI / 3), 0, 0.5;
        stance_foot_endT.p << 0, 0, 0;
    }
    else {
        if (modelPtr_->existJointName("right_toe_roll")) {
            contact_joint_id = modelPtr_->getJointId("right_toe_roll");
        }
        else {
            throw std::runtime_error("Can not find joint: right_toe_roll");
        }

        stance_foot_endT.R << 0,             -1, 0,
                              0.5,           0,  -sin(M_PI / 3),
                              sin(M_PI / 3), 0,  0.5;
        stance_foot_endT.p << 0, 0, 0;
    }

    stance_foot_T_des = stance_foot_T_des_input;
}

void DigitModifiedDynamicsConstraints::reinitialize() {
    // swap the stance leg
    if (stanceLeg == 'L' || stanceLeg == 'l') {
        stanceLeg = 'R';
    }
    else {
        stanceLeg = 'L';
    }

    // reinitialize the stance leg end effector transformation matrix
    if (stanceLeg == 'L' || stanceLeg == 'l') {
        if (modelPtr_->existJointName("left_toe_roll")) {
            contact_joint_id = modelPtr_->getJointId("left_toe_roll");
        }
        else {
            throw std::runtime_error("Can not find joint: left_toe_roll");
        }

        stance_foot_endT.R << 0,             1, 0,
                              -0.5,          0, sin(M_PI / 3),
                              sin(M_PI / 3), 0, 0.5;
        stance_foot_endT.p << 0, 0, 0;
    }
    else {
        if (modelPtr_->existJointName("right_toe_roll")) {
            contact_joint_id = modelPtr_->getJointId("right_toe_roll");
        }
        else {
            throw std::runtime_error("Can not find joint: right_toe_roll");
        }

        stance_foot_endT.R << 0,            -1, 0,
                              0.5,           0, -sin(M_PI / 3),
                              sin(M_PI / 3), 0, 0.5;
        stance_foot_endT.p << 0, 0, 0;
    }
}

int DigitModifiedDynamicsConstraints::return_dependent_joint_index(const int id) {
    assert(0 <= id && id < NUM_DEPENDENT_JOINTS);
    return dependentJointIds[id];
}

int DigitModifiedDynamicsConstraints::return_independent_joint_index(const int id) {
    assert(0 <= id && id < NUM_INDEPENDENT_JOINTS);
    return independentJointIds[id];
}

void DigitModifiedDynamicsConstraints::fill_dependent_vector(VecX& r, const VecX& v, const bool setZero) {
    assert(r.size() == modelPtr_->nv);
    assert(v.size() == NUM_DEPENDENT_JOINTS);

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r(dependentJointIds[i]) = v(i);
    }
}

void DigitModifiedDynamicsConstraints::fill_independent_vector(VecX& r, const VecX& v, const bool setZero) {
    assert(r.size() == modelPtr_->nv);
    assert(v.size() == NUM_INDEPENDENT_JOINTS);

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r(independentJointIds[i]) = v(i);
    }
}

void DigitModifiedDynamicsConstraints::fill_dependent_columns(MatX& r, const MatX& m, const bool setZero) {
    assert(m.cols() == NUM_DEPENDENT_JOINTS);
    assert(r.cols() == modelPtr_->nv);
    assert(m.rows() == r.rows());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.col(dependentJointIds[i]) = m.col(i);
    }
}

void DigitModifiedDynamicsConstraints::fill_independent_columns(MatX& r, const MatX& m, const bool setZero) {
    assert(m.cols() == NUM_INDEPENDENT_JOINTS);
    assert(r.cols() == modelPtr_->nv);
    assert(m.rows() == r.rows());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.col(independentJointIds[i]) = m.col(i);
    }
}

void DigitModifiedDynamicsConstraints::fill_dependent_rows(MatX& r, const MatX& m, const bool setZero) {
    assert(m.rows() == NUM_DEPENDENT_JOINTS);
    assert(r.rows() == modelPtr_->nv);
    assert(m.cols() == r.cols());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.row(dependentJointIds[i]) = m.row(i);
    }
}

void DigitModifiedDynamicsConstraints::fill_independent_rows(MatX& r, const MatX& m, const bool setZero) {
    assert(m.rows() == NUM_INDEPENDENT_JOINTS);
    assert(r.rows() == modelPtr_->nv);
    assert(m.cols() == r.cols());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.row(independentJointIds[i]) = m.row(i);
    }
}

Eigen::VectorXd DigitModifiedDynamicsConstraints::get_dependent_vector(const VecX& v) {
    assert(v.size() == modelPtr_->nv);

    VecX r(NUM_DEPENDENT_JOINTS);

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r(i) = v(dependentJointIds[i]);
    }

    return r;
}

Eigen::VectorXd DigitModifiedDynamicsConstraints::get_independent_vector(const VecX& v) {
    assert(v.size() == modelPtr_->nv);

    VecX r(NUM_INDEPENDENT_JOINTS);

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r(i) = v(independentJointIds[i]);
    }

    return r;
}

void DigitModifiedDynamicsConstraints::get_dependent_columns(MatX& r, const MatX& m) {
    assert(m.cols() == modelPtr_->nv);
    assert(r.cols() == NUM_DEPENDENT_JOINTS);
    assert(m.rows() == r.rows());

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.col(i) = m.col(dependentJointIds[i]);
    }
}

void DigitModifiedDynamicsConstraints::get_independent_columns(MatX& r, const MatX& m) {
    assert(m.cols() == modelPtr_->nv);
    assert(r.cols() == NUM_INDEPENDENT_JOINTS);
    assert(m.rows() == r.rows());

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.col(i) = m.col(independentJointIds[i]);
    }
}

void DigitModifiedDynamicsConstraints::get_dependent_rows(MatX& r, const MatX& m) {
    assert(m.rows() == modelPtr_->nv);
    assert(r.rows() == NUM_DEPENDENT_JOINTS);
    assert(m.cols() == r.cols());

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.row(i) = m.row(dependentJointIds[i]);
    }
}

void DigitModifiedDynamicsConstraints::get_independent_rows(MatX& r, const MatX& m) {
    assert(m.rows() == modelPtr_->nv);
    assert(r.rows() == NUM_INDEPENDENT_JOINTS);
    assert(m.cols() == r.cols());

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.row(i) = m.row(independentJointIds[i]);
    }
}

void DigitModifiedDynamicsConstraints::setupJointPosition(VecX& q, bool compute_derivatives) {
    // fill in dependent joint positions 
    fkhofPtr_->fk(stance_foot_T, *modelPtr_, jtype, contact_joint_id, modelPtr_->getJointId("Rz"), q, stance_foot_endT, startT);
    Transform torso_T = stance_foot_T_des * stance_foot_T.inverse();
    q.block(0, 0, 6, 1) = fkhofPtr_->Transform2xyzrpy(torso_T);

    if (compute_derivatives) {
        get_c(q);
        get_J(q);

        get_dependent_columns(J_dep, J);
        get_independent_columns(J_indep, J);

        J_dep_qr = QRSolver(J_dep);
        J_dep_T_qr = QRSolver(J_dep.transpose());

        // sanity check on uniqueness (these two arguments are actually equivalent)
        if (J_dep_qr.rank() != J_dep.rows() || 
            J_dep_qr.rank() != J_dep.cols() ||
            J_dep_T_qr.rank() != J_dep.rows() || 
            J_dep_T_qr.rank() != J_dep.cols()) {
            throw std::runtime_error("Constraint jacobian is not full rank!");
        }

        P_dep = -J_dep_qr.solve(J_indep);

        pq_dep_pq_indep = P_dep;
    }
}

void DigitModifiedDynamicsConstraints::get_c(const VecX& q) {
    fkhofPtr_->fk(stance_foot_T, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    c = fkhofPtr_->Transform2xyzrpy(stance_foot_T) - fkhofPtr_->Transform2xyzrpy(stance_foot_T_des);
}

void DigitModifiedDynamicsConstraints::get_J(const VecX& q) {
    assert(J.rows() == NUM_DEPENDENT_JOINTS);
    assert(J.cols() == modelPtr_->nv);

    // fkhofPtr_->fk(stance_foot_T, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    fkhofPtr_->fk_jacobian(dTdq, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    fkhofPtr_->Transform2xyzrpyJacobian(J, stance_foot_T, dTdq);
}

void DigitModifiedDynamicsConstraints::get_Jx_partial_dq(const VecX& q, const VecX& x) {
    assert(x.size() == modelPtr_->nv);
    assert(Jx_partial_dq.rows() == NUM_DEPENDENT_JOINTS);
    assert(Jx_partial_dq.cols() == modelPtr_->nv);

    // fkhofPtr_->fk(stance_foot_T, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    // fkhofPtr_->fk_jacobian(dTdq, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    fkhofPtr_->fk_hessian(ddTddq, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    
    Eigen::Array<MatX, 6, 1> H_contact;
    for (int i = 0; i < 6; i++) {
        H_contact(i) = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    }
    fkhofPtr_->Transform2xyzrpyHessian(H_contact, stance_foot_T, dTdq, ddTddq);

    for (int i = 0; i < 6; i++) {
        Jx_partial_dq.row(i) = (H_contact(i) * x).transpose();
    }
}

void DigitModifiedDynamicsConstraints::get_JTx_partial_dq(const VecX& q, const VecX& x) {
    assert(x.size() == NUM_DEPENDENT_JOINTS);
    assert(JTx_partial_dq.rows() == modelPtr_->nv);
    assert(JTx_partial_dq.cols() == modelPtr_->nv);

    // fkhofPtr_->fk(stance_foot_T, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    // fkhofPtr_->fk_jacobian(dTdq, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    fkhofPtr_->fk_hessian(ddTddq, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    
    Eigen::Array<MatX, 6, 1> H_contact;
    for (int i = 0; i < 6; i++) {
        H_contact(i) = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    }
    fkhofPtr_->Transform2xyzrpyHessian(H_contact, stance_foot_T, dTdq, ddTddq);

    JTx_partial_dq.setZero();
    for (int i = 0; i < 6; i++) {
        JTx_partial_dq += H_contact(i) * x(i);
    }
}

void DigitModifiedDynamicsConstraints::get_Jxy_partial_dq(const VecX& q, const VecX& x, const VecX& y) {
    assert(x.size() == modelPtr_->nv);
    assert(y.size() == modelPtr_->nv);
    assert(Jxy_partial_dq.rows() == NUM_DEPENDENT_JOINTS);
    assert(Jxy_partial_dq.cols() == modelPtr_->nv);

    // fkhofPtr_->fk(stance_foot_T, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    // fkhofPtr_->fk_jacobian(dTdq, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    // fkhofPtr_->fk_hessian(ddTddq, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    fkhofPtr_->fk_thirdorder(dddTdddq, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    
    Eigen::Array<MatX, 6, 1> TOx_contact;
    for (int i = 0; i < 6; i++) {
        TOx_contact(i) = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    }
    fkhofPtr_->Transform2xyzrpyThirdOrder(TOx_contact, x, stance_foot_T, dTdq, ddTddq, dddTdddq);

    for (int i = 0; i < 6; i++) {
        Jxy_partial_dq.row(i) = TOx_contact(i) * y;
    }
}

}; // namespace DigitModified
}; // namespace IDTO
