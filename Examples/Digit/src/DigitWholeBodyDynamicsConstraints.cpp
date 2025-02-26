#include "DigitWholeBodyDynamicsConstraints.h"

namespace RAPTOR {
namespace DigitWholeBodySysID {

DigitWholeBodyDynamicsConstraints::DigitWholeBodyDynamicsConstraints(const std::shared_ptr<Model>& modelPtr_input) :
    modelPtr_(modelPtr_input),
    DynamicsConstraints(modelPtr_input->nv, NUM_DEPENDENT_JOINTS + 6) {
    dataPtr_ = std::make_shared<Data>(*modelPtr_);

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

    if (modelPtr_->existFrame("left_foot")) {
        left_foot_idx = modelPtr_->getFrameId("left_foot");
    }
    else {
        throw std::invalid_argument("left foot frame not found!");
    }
    
    if (modelPtr_->existFrame("right_foot")) {
        right_foot_idx = modelPtr_->getFrameId("right_foot");
    }
    else {
        throw std::invalid_argument("right foot frame not found!");
    }

    J_left_foot = MatX::Zero(6, modelPtr_->nv);
    J_right_foot = MatX::Zero(6, modelPtr_->nv);
}

int DigitWholeBodyDynamicsConstraints::return_dependent_joint_index(const int id) {
    assert(0 <= id && id < NUM_DEPENDENT_JOINTS);
    return dependentJointIds[id];
}

int DigitWholeBodyDynamicsConstraints::return_independent_joint_index(const int id) {
    assert(0 <= id && id < NUM_INDEPENDENT_JOINTS);
    return independentJointIds[id];
}

void DigitWholeBodyDynamicsConstraints::fill_dependent_vector(VecX& r, const VecX& v, const bool setZero) {
    assert(r.size() == modelPtr_->nv);
    assert(v.size() == NUM_DEPENDENT_JOINTS);

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r(dependentJointIds[i]) = v(i);
    }
}

void DigitWholeBodyDynamicsConstraints::fill_independent_vector(VecX& r, const VecX& v, const bool setZero) {
    assert(r.size() == modelPtr_->nv);
    assert(v.size() == NUM_INDEPENDENT_JOINTS);

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r(independentJointIds[i]) = v(i);
    }
}

void DigitWholeBodyDynamicsConstraints::fill_dependent_columns(MatX& r, const MatX& m, const bool setZero) {
    assert(m.cols() == NUM_DEPENDENT_JOINTS);
    assert(r.cols() == modelPtr_->nv);
    assert(m.rows() == r.rows());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.col(dependentJointIds[i]) = m.col(i);
    }
}

void DigitWholeBodyDynamicsConstraints::fill_independent_columns(MatX& r, const MatX& m, const bool setZero) {
    assert(m.cols() == NUM_INDEPENDENT_JOINTS);
    assert(r.cols() == modelPtr_->nv);
    assert(m.rows() == r.rows());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.col(independentJointIds[i]) = m.col(i);
    }
}

void DigitWholeBodyDynamicsConstraints::fill_dependent_rows(MatX& r, const MatX& m, const bool setZero) {
    assert(m.rows() == NUM_DEPENDENT_JOINTS);
    assert(r.rows() == modelPtr_->nv);
    assert(m.cols() == r.cols());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.row(dependentJointIds[i]) = m.row(i);
    }
}

void DigitWholeBodyDynamicsConstraints::fill_independent_rows(MatX& r, const MatX& m, const bool setZero) {
    assert(m.rows() == NUM_INDEPENDENT_JOINTS);
    assert(r.rows() == modelPtr_->nv);
    assert(m.cols() == r.cols());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.row(independentJointIds[i]) = m.row(i);
    }
}

Eigen::VectorXd DigitWholeBodyDynamicsConstraints::get_dependent_vector(const VecX& v) {
    assert(v.size() == modelPtr_->nv);

    VecX r(NUM_DEPENDENT_JOINTS);

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r(i) = v(dependentJointIds[i]);
    }

    return r;
}

Eigen::VectorXd DigitWholeBodyDynamicsConstraints::get_independent_vector(const VecX& v) {
    assert(v.size() == modelPtr_->nv);

    VecX r(NUM_INDEPENDENT_JOINTS);

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r(i) = v(independentJointIds[i]);
    }

    return r;
}

void DigitWholeBodyDynamicsConstraints::get_dependent_columns(MatX& r, const MatX& m) {
    assert(m.cols() == modelPtr_->nv);
    
    r.resize(m.rows(), NUM_DEPENDENT_JOINTS);

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.col(i) = m.col(dependentJointIds[i]);
    }
}

void DigitWholeBodyDynamicsConstraints::get_independent_columns(MatX& r, const MatX& m) {
    assert(m.cols() == modelPtr_->nv);
    
    r.resize(m.rows(), NUM_INDEPENDENT_JOINTS);

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.col(i) = m.col(independentJointIds[i]);
    }
}

void DigitWholeBodyDynamicsConstraints::get_dependent_rows(MatX& r, const MatX& m) {
    assert(m.rows() == modelPtr_->nv);

    r.resize(NUM_DEPENDENT_JOINTS, m.cols());

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.row(i) = m.row(dependentJointIds[i]);
    }
}

void DigitWholeBodyDynamicsConstraints::get_independent_rows(MatX& r, const MatX& m) {
    assert(m.rows() == modelPtr_->nv);

    r.resize(NUM_INDEPENDENT_JOINTS, m.cols());

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.row(i) = m.row(independentJointIds[i]);
    }
}

void DigitWholeBodyDynamicsConstraints::get_J(const VecX& q) {
    assert(J.rows() == NUM_DEPENDENT_JOINTS + 6);
    assert(J.cols() == modelPtr_->nv);

    J.setZero();

    const double& q1 = q(0);
    const double& q2 = q(1);
    const double& q3 = q(2);
    const double& q4 = q(3);
    const double& q5 = q(4);
    const double& q6 = q(5);
    const double& q7 = q(6);
    const double& q8 = q(7);
    const double& q9 = q(8);
    const double& q10 = q(9);
    const double& q11 = q(10);
    const double& q12 = q(11);
    const double& q13 = q(12);
    const double& q14 = q(13);
    const double& q15 = q(14);
    const double& q16 = q(15);
    const double& q17 = q(16);
    const double& q18 = q(17);
    const double& q19 = q(18);
    const double& q20 = q(19);
    const double& q21 = q(20);
    const double& q22 = q(21);
    const double& q23 = q(22);
    const double& q24 = q(23);
    const double& q25 = q(24);
    const double& q26 = q(25);
    const double& q27 = q(26);
    const double& q28 = q(27);
    const double& q29 = q(28);
    const double& q30 = q(29);
    const double& q31 = q(30);
    const double& q32 = q(31);
    const double& q33 = q(32);
    const double& q34 = q(33);
    const double& q35 = q(34);
    const double& q36 = q(35);
    const double& q37 = q(36);
    const double& q38 = q(37);
    const double& q39 = q(38);
    const double& q40 = q(39);
    const double& q41 = q(40);
    const double& q42 = q(41);
    const double& q43 = q(42);
    const double& q44 = q(43);

    double t2 = cos(q10);
    double t3 = cos(q11);
    double t4 = sin(q10);
    double t5 = sin(q11);
    double t6 = q12+q13-5.308372276804118E-4;
    double t7 = q12+q13+1.570265489567216;
    double t8 = cos(t7);
    double t9 = cos(t6);
    double t10 = t9*1.200056415552922E-1;
    double t11 = t8*1.200056415552922E-1;
    J(0, 9) = t3*cos(q10+1.349139511789001E-1)*4.999998096141249E-1;
    J(0, 10) = t3*-4.363322573663341E-4-t2*t5*6.725249670103779E-2-t4*t5*4.95456265780985E-1;
    J(0, 11) = t11-cos(q12+1.352236194373202E-1)*5.000003989318408E-1;
    J(0, 12) = t11;
    J(1, 9) = t3*cos(q10-1.435819975335223)*4.98726624169518E-1;
    J(1, 10) = t3*3.566166491453071E-2+t2*t5*4.941904630377533E-1-t4*t5*6.711208459029357E-2;
    J(1, 11) = t10-cos(q12-1.435572707357576)*5.000003989318408E-1;
    J(1, 12) = t10;
    J(2, 9) = t3*cos(q10-1.448085935274989)*3.566433414933032E-2;
    J(2, 10) = t3*-4.987264332975392E-1+t2*t5*3.5396156994904E-2-t4*t5*4.36540952342363E-3;

    t2 = cos(q14);
    t3 = cos(q15);
    t4 = cos(q16);
    t5 = cos(q17);
    t6 = cos(q18);
    t7 = cos(q19);
    t8 = sin(q14);
    t9 = sin(q15);
    t10 = sin(q16);
    t11 = sin(q17);
    double t12 = sin(q18);
    double t13 = sin(q19);
    double t14 = q20+1.252274463839301;
    double t15 = q20-3.185218629555953E-1;
    double t18 = q20+q21+1.670234059733348E-1;
    double t19 = q20+q21-1.403772920821562;
    double t20 = q20+q21+7.621286815357554E-1;
    double t21 = q20+q21-8.086676452591412E-1;
    double t22 = -q20;
    double t23 = -q21;
    double t16 = cos(t14);
    double t17 = cos(t15);
    double t24 = cos(t20);
    double t25 = cos(t21);
    double t26 = cos(t18);
    double t27 = cos(t19);
    double t31 = q20+t23+7.667291949103712E-1;
    double t32 = q21+t22+8.040671318845255E-1;
    double t40 = q20+t23+1.716239193479506E-1;
    double t41 = q21+t22+1.399172407446946;
    double t28 = t16*5.416399999999999E-2;
    double t29 = t17*-5.416399999999999E-2;
    double t30 = t17*5.416399999999999E-2;
    double t33 = t26*-1.023268783116147E-2;
    double t34 = t26*1.023268783116147E-2;
    double t35 = t27*-1.023268783116147E-2;
    double t36 = t27*1.023268783116147E-2;
    double t37 = t24*-1.014435312131828E-2;
    double t38 = t24*1.014435312131828E-2;
    double t39 = t25*1.014435312131828E-2;
    double t42 = cos(t40);
    double t43 = cos(t41);
    double t44 = cos(t31);
    double t45 = cos(t32);
    double t46 = t42*1.014435312131828E-2;
    double t47 = t43*1.014435312131828E-2;
    double t48 = t44*1.023268783116147E-2;
    double t49 = t45*1.023268783116147E-2;
    J(3, 13) = t2*-5.69906102430749E-2-t8*1.034574367518647E-3+t2*t10*1.265737725314828E-2+t8*t10*2.246288488292063E-3-t2*t3*t4*3.020093502477822E-3-t2*t4*t9*3.397508937976584E-1-t3*t4*t8*3.399783937080488E-1+t4*t8*t9*3.105800979186457E-3;
    J(3, 14) = t4*(cos(q14+q15-1.561784314582967)*3.398784455052338E+34+cos(q14-q15+1.210504281999349)*1.215544947661321E+31)*(-1.0E-35);
    J(3, 15) = t2*t4*-2.246288488292063E-3+t4*t8*1.265737725314828E-2-t2*t3*t10*3.399783937080488E-1+t2*t9*t10*3.105800979186457E-3+t3*t8*t10*3.020093502477822E-3+t8*t9*t10*3.397508937976584E-1;
    J(3, 19) = t29+t39-t46;
    J(3, 20) = t39+t46;
    J(4, 13) = t2*-1.034574367518647E-3+t8*5.69906102430749E-2+t2*t10*2.246288488292063E-3-t8*t10*1.265737725314828E-2-t2*t3*t4*3.399783937080488E-1+t2*t4*t9*3.105800979186457E-3+t3*t4*t8*3.020093502477822E-3+t4*t8*t9*3.397508937976584E-1;
    J(4, 14) = t4*(cos(-q14+q15+3.602920447955478E-1)*1.215544947661321E+31-cos(q14+q15+9.0120122119293E-3)*3.398784455052338E+34)*1.0E-35;
    J(4, 15) = t2*t4*1.265737725314828E-2+t4*t8*2.246288488292063E-3+t2*t3*t10*3.020093502477822E-3+t2*t9*t10*3.397508937976584E-1+t3*t8*t10*3.399783937080488E-1-t8*t9*t10*3.105800979186457E-3;
    J(4, 19) = t28+t37-t47;
    J(4, 20) = t37+t47;
    J(5, 14) = t4*cos(q15+1.846520285037398E-1)*1.285515503217091E-2;
    J(5, 15) = t4*3.397568910104677E-1-t3*t10*2.360264165562881E-3-t9*t10*1.263661995827637E-2;
    J(5, 20) = cos(q21+1.080650544491351)*2.028870624263657E-2;
    J(6, 16) = t5*5.699060997402856E-2-t11*1.034589188110661E-3+t5*t13*1.549180159908108E-2+t11*t13*8.23479231089269E-4+t5*t6*t7*3.104080344633564E-3-t5*t7*t12*2.87566285868891E-1-t6*t7*t11*2.879825211612492E-1-t7*t11*t12*3.064210757541306E-3;
    J(6, 17) = t7*(cos(q17+q18+1.56007950333587)*1.438954648974688E+34-cos(q17-q18+1.475301494119507)*1.045351025312074E+31)*2.0E-35;
    J(6, 18) = t5*t7*-8.23479231089269E-4+t7*t11*1.549180159908108E-2-t5*t6*t13*2.879825211612492E-1-t5*t12*t13*3.064210757541306E-3-t6*t11*t13*3.104080344633564E-3+t11*t12*t13*2.87566285868891E-1;
    J(6, 19) = t29+t33+t49;
    J(6, 20) = t33-t49;
    J(7, 16) = t5*1.034589188110661E-3+t11*5.699060997402856E-2-t5*t13*8.23479231089269E-4+t11*t13*1.549180159908108E-2+t5*t6*t7*2.879825211612492E-1+t5*t7*t12*3.064210757541306E-3+t6*t7*t11*3.104080344633564E-3-t7*t11*t12*2.87566285868891E-1;
    J(7, 17) = t7*(cos(-q17+q18+9.549483267538958E-2)*1.045351025312074E+31-cos(q17+q18-1.071682345902628E-2)*1.438954648974688E+34)*(-2.0E-35);
    J(7, 18) = t5*t7*-1.549180159908108E-2-t7*t11*8.23479231089269E-4+t5*t6*t13*3.104080344633564E-3-t5*t12*t13*2.87566285868891E-1-t6*t11*t13*2.879825211612492E-1-t11*t12*t13*3.064210757541306E-3;
    J(7, 19) = t28+t35-t48;
    J(7, 20) = t35+t48;
    J(8, 17) = t7*cos(q18+4.2389004608181E-2)*-1.551367251263628E-2;
    J(8, 18) = t7*-2.875818595898751E-1+t6*t13*6.574122182670536E-4+t12*t13*1.549973690114125E-2;
    J(8, 20) = cos(q21-1.085251057865967)*2.046537566232294E-2;

    t2 = cos(q29);
    t3 = cos(q30);
    t4 = sin(q29);
    t5 = sin(q30);
    t6 = q31+q32+5.308372276804118E-4;
    t7 = q31+q32-1.570265489567216;
    t8 = cos(t7);
    t9 = cos(t6);
    t10 = t9*1.200056415552922E-1;
    t11 = t8*-1.200056415552922E-1;
    t12 = t8*1.200056415552922E-1;
    J(9, 28) = t3*cos(q29-1.349132376885636E-1)*-4.999998093925668E-1;
    J(9, 29) = t3*-4.365860704565479E-4-t2*t5*6.725214316796235E-2+t4*t5*4.954563135453205E-1;
    J(9, 30) = t11+cos(q31-1.352236194373202E-1)*5.000003989318408E-1;
    J(9, 31) = t11;
    J(10, 28) = t3*cos(q29+1.435820652758106)*4.987266335404373E-1;
    J(10, 29) = t3*-3.566153386244494E-2-t2*t5*4.941905177865888E-1-t4*t5*6.711175107535901E-2;
    J(10, 30) = t10-cos(q31+1.435572707357576)*5.000003989318408E-1;
    J(10, 31) = t10;
    J(11, 28) = t3*cos(q29+1.448093791924176)*-3.566420621322202E-2;
    J(11, 29) = t3*-4.987264424463383E-1+t2*t5*3.539606431708407E-2+t4*t5*4.365115769377903E-3;

    t2 = cos(q33);
    t3 = cos(q34);
    t4 = cos(q35);
    t5 = cos(q36);
    t6 = cos(q37);
    t7 = cos(q38);
    t8 = sin(q33);
    t9 = sin(q34);
    t10 = sin(q35);
    t11 = sin(q36);
    t12 = sin(q37);
    t13 = sin(q38);
    t14 = q39-1.252274463839301;
    t15 = q39+3.185218629555953E-1;
    t18 = q39+q40-1.670234059733348E-1;
    t19 = q39+q40+1.403772920821562;
    t20 = q39+q40-7.621286815357554E-1;
    t21 = q39+q40+8.086676452591412E-1;
    t22 = -q39;
    t23 = -q40;
    t16 = cos(t14);
    t17 = cos(t15);
    t24 = cos(t20);
    t25 = cos(t21);
    t26 = cos(t18);
    t27 = cos(t19);
    t30 = q40+t22+7.667291949103712E-1;
    t31 = q39+t23+8.040671318845255E-1;
    t39 = q40+t22+1.716239193479506E-1;
    t40 = q39+t23+1.399172407446946;
    t28 = t16*5.416399999999999E-2;
    t29 = t17*5.416399999999999E-2;
    t32 = t26*1.023268783116147E-2;
    t33 = t27*-1.023268783116147E-2;
    t34 = t27*1.023268783116147E-2;
    t35 = t24*-1.014435312131828E-2;
    t36 = t24*1.014435312131828E-2;
    t37 = t25*-1.014435312131828E-2;
    t38 = t25*1.014435312131828E-2;
    t41 = cos(t39);
    t42 = cos(t40);
    t43 = cos(t30);
    t44 = cos(t31);
    t45 = t41*1.014435312131828E-2;
    t46 = t42*1.014435312131828E-2;
    t47 = t43*1.023268783116147E-2;
    t48 = t44*1.023268783116147E-2;
    J(12, 32) = t2*5.699060997402857E-2-t8*1.034589188110661E-3-t2*t10*1.26575449899492E-2+t8*t10*2.246221860400801E-3+t2*t3*t4*3.020283789547076E-3-t2*t4*t9*3.397508858570615E-1-t3*t4*t8*3.399783924207052E-1-t4*t8*t9*3.105990081579733E-3;
    J(12, 33) = t4*(cos(q33+q34+1.56178375635807)*3.39878442601012E+34+cos(-q33+q34+1.210518490104648)*1.215573989880305E+31)*1.0E-35;
    J(12, 34) = t2*t4*-2.246221860400801E-3-t4*t8*1.26575449899492E-2-t2*t3*t10*3.399783924207052E-1-t2*t9*t10*3.105990081579733E-3-t3*t8*t10*3.020283789547076E-3+t8*t9*t10*3.397508858570615E-1;
    J(12, 38) = t29+t37+t45;
    J(12, 39) = t37-t45;
    J(13, 32) = t2*-1.034589188110661E-3-t8*5.699060997402857E-2+t2*t10*2.246221860400801E-3+t8*t10*1.26575449899492E-2-t2*t3*t4*3.399783924207052E-1-t2*t4*t9*3.105990081579733E-3-t3*t4*t8*3.020283789547076E-3+t4*t8*t9*3.397508858570615E-1;
    J(13, 33) = (t4*(cos(q33-q34+3.602778366902485E-1)*1.215573989880305E+31-cos(q33+q34-9.01257043682671E-3)*3.39878442601012E+34))/1.0E+35;
    J(13, 34) = t2*t4*-1.26575449899492E-2+t4*t8*2.246221860400801E-3-t2*t3*t10*3.020283789547076E-3+t2*t9*t10*3.397508858570615E-1+t3*t8*t10*3.399783924207052E-1+t8*t9*t10*3.105990081579733E-3;
    J(13, 38) = t28+t35-t46;
    J(13, 39) = t35+t46;
    J(14, 33) = t4*cos(q34-1.846452035635318E-1)*-1.28553085462283E-2;
    J(14, 34) = t4*3.39756885202024E-1-t3*t10*2.360206106172353E-3+t9*t10*1.263678697118524E-2;
    J(14, 39) = cos(q40-1.080650544491351)*-2.028870624263657E-2;
    J(15, 35) = t5*-5.699060997402856E-2-t11*1.034589188110661E-3-t5*t13*1.549180159908108E-2+t11*t13*8.23479231089269E-4-t5*t6*t7*3.104080344633564E-3-t5*t7*t12*2.87566285868891E-1-t6*t7*t11*2.879825211612492E-1+t7*t11*t12*3.064210757541306E-3;
    J(15, 36) = t7*(cos(q36+q37-1.56007950333587)*1.438954648974688E+34-cos(-q36+q37+1.475301494119507)*1.045351025312074E+31)*(-2.0E-35);
    J(15, 37) = t5*t7*-8.23479231089269E-4-t7*t11*1.549180159908108E-2-t5*t6*t13*2.879825211612492E-1+t5*t12*t13*3.064210757541306E-3+t6*t11*t13*3.104080344633564E-3+t11*t12*t13*2.87566285868891E-1;
    J(15, 38) = t29+t32-t48;
    J(15, 39) = t32+t48;
    J(16, 35) = t5*1.034589188110661E-3-t11*5.699060997402856E-2-t5*t13*8.23479231089269E-4-t11*t13*1.549180159908108E-2+t5*t6*t7*2.879825211612492E-1-t5*t7*t12*3.064210757541306E-3-t6*t7*t11*3.104080344633564E-3-t7*t11*t12*2.87566285868891E-1;
    J(16, 36) = t7*(cos(q36-q37+9.549483267538958E-2)*1.045351025312074E+31-cos(q36+q37+1.071682345902628E-2)*1.438954648974688E+34)*-2.0E-35;
    J(16, 37) = t5*t7*1.549180159908108E-2-t7*t11*8.23479231089269E-4-t5*t6*t13*3.104080344633564E-3-t5*t12*t13*2.87566285868891E-1-t6*t11*t13*2.879825211612492E-1+t11*t12*t13*3.064210757541306E-3;
    J(16, 38) = t28+t33-t47;
    J(16, 39) = t33+t47;
    J(17, 36) = t7*cos(q37-4.2389004608181E-2)*1.551367251263628E-2;
    J(17, 37) = t7*-2.875818595898751E-1+t6*t13*6.574122182670536E-4-t12*t13*1.549973690114125E-2;
    J(17, 39) = cos(q40+1.085251057865967)*-2.046537566232294E-2;

    pinocchio::forwardKinematics(*modelPtr_, *dataPtr_, q);
    pinocchio::computeJointJacobians(*modelPtr_, *dataPtr_, q);
    pinocchio::updateFramePlacements(*modelPtr_, *dataPtr_);

    // left foot jacobian
    pinocchio::getFrameJacobian(
        *modelPtr_, *dataPtr_, 
        left_foot_idx, 
        pinocchio::LOCAL_WORLD_ALIGNED,
        J_left_foot);
    J.middleRows(18, 6) = J_left_foot;

    // right foot jacobian
    pinocchio::getFrameJacobian(
        *modelPtr_, *dataPtr_, 
        right_foot_idx, 
        pinocchio::LOCAL_WORLD_ALIGNED,
        J_right_foot);
    J.middleRows(24, 6) = J_right_foot;
}

}; // namespace DigitWholeBodySysID
}; // namespace RAPTOR
