#include "TalosDynamicsConstraints.h"

namespace RAPTOR {
namespace Talos {

TalosDynamicsConstraints::TalosDynamicsConstraints(const std::shared_ptr<Model>& modelPtr_input, 
                                                   char stanceLeg_input, 
                                                   const Transform& stance_foot_T_des_input) :
    modelPtr_(modelPtr_input),
    stanceLeg(stanceLeg_input),
    DynamicsConstraints(modelPtr_input->nv, NUM_DEPENDENT_JOINTS) {
    fkPtr_ = std::make_unique<ForwardKinematicsSolver>(modelPtr_.get());

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
            throw std::runtime_error("Can not find joint: " + independentJointNames[i]);
        }
    }

    if (stanceLeg == 'L' || stanceLeg == 'l') {
        contact_joint_id = modelPtr_->getJointId(std::string(LEFT_FOOT_NAME));
        stance_foot_endT.p << 0, 0, -0.107;
    }
    else {
        contact_joint_id = modelPtr_->getJointId(std::string(RIGHT_FOOT_NAME));
        stance_foot_endT.p << 0, 0, -0.107;
    }

    stance_foot_T_des = stance_foot_T_des_input;
}

void TalosDynamicsConstraints::reinitialize() {
    // swap the stance leg
    if (stanceLeg == 'L' || stanceLeg == 'l') {
        stanceLeg = 'R';
    }
    else {
        stanceLeg = 'L';
    }

    // reinitialize the stance leg end effector transformation matrix
    if (stanceLeg == 'L' || stanceLeg == 'l') {
        contact_joint_id = modelPtr_->getJointId(std::string(LEFT_FOOT_NAME));
        stance_foot_endT.p << 0, 0, -0.107;
    }
    else {
        contact_joint_id = modelPtr_->getJointId(std::string(RIGHT_FOOT_NAME));
        stance_foot_endT.p << 0, 0, -0.107;
    }
}

int TalosDynamicsConstraints::return_dependent_joint_index(const int id) {
    assert(0 <= id && id < NUM_DEPENDENT_JOINTS);
    return dependentJointIds[id];
}

int TalosDynamicsConstraints::return_independent_joint_index(const int id) {
    assert(0 <= id && id < NUM_INDEPENDENT_JOINTS);
    return independentJointIds[id];
}

void TalosDynamicsConstraints::fill_dependent_vector(VecX& r, const VecX& v, const bool setZero) {
    assert(r.size() == modelPtr_->nv);
    assert(v.size() == NUM_DEPENDENT_JOINTS);

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r(dependentJointIds[i]) = v(i);
    }
}

void TalosDynamicsConstraints::fill_independent_vector(VecX& r, const VecX& v, const bool setZero) {
    assert(r.size() == modelPtr_->nv);
    assert(v.size() == NUM_INDEPENDENT_JOINTS);

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r(independentJointIds[i]) = v(i);
    }
}

void TalosDynamicsConstraints::fill_dependent_columns(MatX& r, const MatX& m, const bool setZero) {
    assert(m.cols() == NUM_DEPENDENT_JOINTS);
    assert(r.cols() == modelPtr_->nv);
    assert(m.rows() == r.rows());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.col(dependentJointIds[i]) = m.col(i);
    }
}

void TalosDynamicsConstraints::fill_independent_columns(MatX& r, const MatX& m, const bool setZero) {
    assert(m.cols() == NUM_INDEPENDENT_JOINTS);
    assert(r.cols() == modelPtr_->nv);
    assert(m.rows() == r.rows());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.col(independentJointIds[i]) = m.col(i);
    }
}

void TalosDynamicsConstraints::fill_dependent_rows(MatX& r, const MatX& m, const bool setZero) {
    assert(m.rows() == NUM_DEPENDENT_JOINTS);
    assert(r.rows() == modelPtr_->nv);
    assert(m.cols() == r.cols());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.row(dependentJointIds[i]) = m.row(i);
    }
}

void TalosDynamicsConstraints::fill_independent_rows(MatX& r, const MatX& m, const bool setZero) {
    assert(m.rows() == NUM_INDEPENDENT_JOINTS);
    assert(r.rows() == modelPtr_->nv);
    assert(m.cols() == r.cols());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.row(independentJointIds[i]) = m.row(i);
    }
}

Eigen::VectorXd TalosDynamicsConstraints::get_dependent_vector(const VecX& v) {
    assert(v.size() == modelPtr_->nv);

    VecX r(NUM_DEPENDENT_JOINTS);

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r(i) = v(dependentJointIds[i]);
    }

    return r;
}

Eigen::VectorXd TalosDynamicsConstraints::get_independent_vector(const VecX& v) {
    assert(v.size() == modelPtr_->nv);

    VecX r(NUM_INDEPENDENT_JOINTS);

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r(i) = v(independentJointIds[i]);
    }

    return r;
}

void TalosDynamicsConstraints::get_dependent_columns(MatX& r, const MatX& m) {
    assert(m.cols() == modelPtr_->nv);
    assert(r.cols() == NUM_DEPENDENT_JOINTS);
    assert(m.rows() == r.rows());

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.col(i) = m.col(dependentJointIds[i]);
    }
}

void TalosDynamicsConstraints::get_independent_columns(MatX& r, const MatX& m) {
    assert(m.cols() == modelPtr_->nv);
    assert(r.cols() == NUM_INDEPENDENT_JOINTS);
    assert(m.rows() == r.rows());

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.col(i) = m.col(independentJointIds[i]);
    }
}

void TalosDynamicsConstraints::get_dependent_rows(MatX& r, const MatX& m) {
    assert(m.rows() == modelPtr_->nv);
    assert(r.rows() == NUM_DEPENDENT_JOINTS);
    assert(m.cols() == r.cols());

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.row(i) = m.row(dependentJointIds[i]);
    }
}

void TalosDynamicsConstraints::get_independent_rows(MatX& r, const MatX& m) {
    assert(m.rows() == modelPtr_->nv);
    assert(r.rows() == NUM_INDEPENDENT_JOINTS);
    assert(m.cols() == r.cols());

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.row(i) = m.row(independentJointIds[i]);
    }
}

void TalosDynamicsConstraints::setupJointPosition(VecX& q, bool compute_derivatives) {
    if (recoverSavedData(q, compute_derivatives)) {
        return;
    }

    // fill in dependent joint positions 
    fkPtr_->compute(modelPtr_->getJointId("Rz"), 
                    contact_joint_id, 
                    q, 
                    nullptr, 
                    &stance_foot_endT, 
                    0);

    Transform torso_T = stance_foot_T_des * fkPtr_->getTransform().inverse();
    q.head(6) = torso_T.getXYZRPY();
    
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

    updateQueue(q, compute_derivatives);      
}

void TalosDynamicsConstraints::get_c(const VecX& q) {
    fkPtr_->compute(0, 
                    contact_joint_id, 
                    q, 
                    nullptr, 
                    &stance_foot_endT, 
                    0);

    c_translation = fkPtr_->getTranslation() - stance_foot_T_des.p;
    c_rotation = fkPtr_->getRPY() - stance_foot_T_des.getRPY();

    c = VecX(6);
    c << c_translation, c_rotation;
}

void TalosDynamicsConstraints::get_J(const VecX& q) {
    assert(J.rows() == NUM_DEPENDENT_JOINTS);
    assert(J.cols() == modelPtr_->nv);

    fkPtr_->compute(0, 
                    contact_joint_id, 
                    q, 
                    nullptr, 
                    &stance_foot_endT, 
                    1);

    J_translation = fkPtr_->getTranslationJacobian();
    J_rotation = fkPtr_->getRPYJacobian();

    J << J_translation, J_rotation;
}

void TalosDynamicsConstraints::get_Jx_partial_dq(const VecX& q, const VecX& x) {
    assert(x.size() == modelPtr_->nv);
    assert(Jx_partial_dq.rows() == NUM_DEPENDENT_JOINTS);
    assert(Jx_partial_dq.cols() == modelPtr_->nv);

    fkPtr_->compute(0, 
                    contact_joint_id, 
                    q, 
                    nullptr, 
                    &stance_foot_endT, 
                    2);

    fkPtr_->getTranslationHessian(H_translation);
    fkPtr_->getRPYHessian(H_rotation);

    Jx_partial_dq.row(0) = H_translation(0) * x;
    Jx_partial_dq.row(1) = H_translation(1) * x;
    Jx_partial_dq.row(2) = H_translation(2) * x;
    Jx_partial_dq.row(3) = H_rotation(0) * x;
    Jx_partial_dq.row(4) = H_rotation(1) * x;
    Jx_partial_dq.row(5) = H_rotation(2) * x;
}

void TalosDynamicsConstraints::get_JTx_partial_dq(const VecX& q, const VecX& x) {
    assert(x.size() == NUM_DEPENDENT_JOINTS);
    assert(JTx_partial_dq.rows() == modelPtr_->nv);
    assert(JTx_partial_dq.cols() == modelPtr_->nv);

    fkPtr_->compute(0, 
                    contact_joint_id, 
                    q, 
                    nullptr, 
                    &stance_foot_endT, 
                    2);

    fkPtr_->getTranslationHessian(H_translation);
    fkPtr_->getRPYHessian(H_rotation);

    JTx_partial_dq.setZero();
    JTx_partial_dq += H_translation(0) * x(0);
    JTx_partial_dq += H_translation(1) * x(1);
    JTx_partial_dq += H_translation(2) * x(2);
    JTx_partial_dq += H_rotation(0) * x(3);
    JTx_partial_dq += H_rotation(1) * x(4);
    JTx_partial_dq += H_rotation(2) * x(5);
}

void TalosDynamicsConstraints::get_Jxy_partial_dq(const VecX& q, const VecX& x, const VecX& y) {
    assert(x.size() == modelPtr_->nv);
    assert(y.size() == modelPtr_->nv);
    assert(Jxy_partial_dq.rows() == NUM_DEPENDENT_JOINTS);
    assert(Jxy_partial_dq.cols() == modelPtr_->nv);

    fkPtr_->compute(0, 
                    contact_joint_id, 
                    q, 
                    nullptr, 
                    &stance_foot_endT, 
                    3);

    fkPtr_->getTranslationThirdOrderTensor(x, TOx_translation);
    fkPtr_->getRPYThirdOrderTensor(x, TOx_rotation);

    Jxy_partial_dq.row(0) = TOx_translation(0) * y;
    Jxy_partial_dq.row(1) = TOx_translation(1) * y;
    Jxy_partial_dq.row(2) = TOx_translation(2) * y;
    Jxy_partial_dq.row(3) = TOx_rotation(0) * y;
    Jxy_partial_dq.row(4) = TOx_rotation(1) * y;
    Jxy_partial_dq.row(5) = TOx_rotation(2) * y;
}

}; // namespace DigitModified
}; // namespace RAPTOR
