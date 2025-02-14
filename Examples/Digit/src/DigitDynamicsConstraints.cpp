#include "DigitDynamicsConstraints.h"

namespace RAPTOR {
namespace Digit {

DigitDynamicsConstraints::DigitDynamicsConstraints(const std::shared_ptr<Model>& modelPtr_input, 
                                                   const char stanceLeg_input, 
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
        if (modelPtr_->existJointName(std::string(LEFT_FOOT_NAME))) {
            contact_joint_id = modelPtr_->getJointId(std::string(LEFT_FOOT_NAME));
        }
        else {
            throw std::runtime_error("Can not find joint: " + std::string(LEFT_FOOT_NAME));
        }

        stance_foot_endT.R << 0,             1, 0,
                              -0.5,          0, sin(M_PI / 3),
                              sin(M_PI / 3), 0, 0.5;
        stance_foot_endT.p << 0, -0.05456, -0.0315;
    }
    else {
        if (modelPtr_->existJointName(std::string(RIGHT_FOOT_NAME))) {
            contact_joint_id = modelPtr_->getJointId(std::string(RIGHT_FOOT_NAME));
        }
        else {
            throw std::runtime_error("Can not find joint: " + std::string(RIGHT_FOOT_NAME));
        }

        stance_foot_endT.R << 0,             -1, 0,
                              0.5,           0,  -sin(M_PI / 3),
                              sin(M_PI / 3), 0,  0.5;
        stance_foot_endT.p << 0, 0.05456, -0.0315;
    }

    stance_foot_T_des = stance_foot_T_des_input;
}

void DigitDynamicsConstraints::reinitialize(const char stanceLeg_input) {
    if (stanceLeg_input == 0) { // swap the stance leg if there's no input
        if (stanceLeg == 'L' || stanceLeg == 'l') {
            stanceLeg = 'R';
        }
        else {
            stanceLeg = 'L';
        }
    }
    else {
        if (stanceLeg_input != 'L' && 
            stanceLeg_input != 'R' && 
            stanceLeg_input != 'l' && 
            stanceLeg_input != 'r') {
            throw std::runtime_error("Invalid stance leg input");
        }
        stanceLeg = stanceLeg_input;
    }

    // reinitialize the stance leg end effector transformation matrix
    if (stanceLeg == 'L' || stanceLeg == 'l') {
        contact_joint_id = modelPtr_->getJointId(std::string(LEFT_FOOT_NAME));
        stance_foot_endT.R << 0,             1, 0,
                              -0.5,          0, sin(M_PI / 3),
                              sin(M_PI / 3), 0, 0.5;
        stance_foot_endT.p << 0, -0.05456, -0.0315;
    }
    else {
        contact_joint_id = modelPtr_->getJointId(std::string(RIGHT_FOOT_NAME));
        stance_foot_endT.R << 0,             -1, 0,
                              0.5,           0,  -sin(M_PI / 3),
                              sin(M_PI / 3), 0,  0.5;
        stance_foot_endT.p << 0, 0.05456, -0.0315;
    }
}

int DigitDynamicsConstraints::return_dependent_joint_index(const int id) {
    assert(0 <= id && id < NUM_DEPENDENT_JOINTS);
    return dependentJointIds[id];
}

int DigitDynamicsConstraints::return_independent_joint_index(const int id) {
    assert(0 <= id && id < NUM_INDEPENDENT_JOINTS);
    return independentJointIds[id];
}

void DigitDynamicsConstraints::fill_dependent_vector(VecX& r, const VecX& v, const bool setZero) {
    assert(r.size() == modelPtr_->nv);
    assert(v.size() == NUM_DEPENDENT_JOINTS);

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r(dependentJointIds[i]) = v(i);
    }
}

void DigitDynamicsConstraints::fill_independent_vector(VecX& r, const VecX& v, const bool setZero) {
    assert(r.size() == modelPtr_->nv);
    assert(v.size() == NUM_INDEPENDENT_JOINTS);

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r(independentJointIds[i]) = v(i);
    }
}

void DigitDynamicsConstraints::fill_dependent_columns(MatX& r, const MatX& m, const bool setZero) {
    assert(m.cols() == NUM_DEPENDENT_JOINTS);
    assert(r.cols() == modelPtr_->nv);
    assert(m.rows() == r.rows());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.col(dependentJointIds[i]) = m.col(i);
    }
}

void DigitDynamicsConstraints::fill_independent_columns(MatX& r, const MatX& m, const bool setZero) {
    assert(m.cols() == NUM_INDEPENDENT_JOINTS);
    assert(r.cols() == modelPtr_->nv);
    assert(m.rows() == r.rows());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.col(independentJointIds[i]) = m.col(i);
    }
}

void DigitDynamicsConstraints::fill_dependent_rows(MatX& r, const MatX& m, const bool setZero) {
    assert(m.rows() == NUM_DEPENDENT_JOINTS);
    assert(r.rows() == modelPtr_->nv);
    assert(m.cols() == r.cols());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.row(dependentJointIds[i]) = m.row(i);
    }
}

void DigitDynamicsConstraints::fill_independent_rows(MatX& r, const MatX& m, const bool setZero) {
    assert(m.rows() == NUM_INDEPENDENT_JOINTS);
    assert(r.rows() == modelPtr_->nv);
    assert(m.cols() == r.cols());

    if (setZero) r.setZero();

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.row(independentJointIds[i]) = m.row(i);
    }
}

Eigen::VectorXd DigitDynamicsConstraints::get_dependent_vector(const VecX& v) {
    assert(v.size() == modelPtr_->nv);

    VecX r(NUM_DEPENDENT_JOINTS);

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r(i) = v(dependentJointIds[i]);
    }

    return r;
}

Eigen::VectorXd DigitDynamicsConstraints::get_independent_vector(const VecX& v) {
    assert(v.size() == modelPtr_->nv);

    VecX r(NUM_INDEPENDENT_JOINTS);

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r(i) = v(independentJointIds[i]);
    }

    return r;
}

void DigitDynamicsConstraints::get_dependent_columns(MatX& r, const MatX& m) {
    assert(m.cols() == modelPtr_->nv);
    
    r.resize(m.rows(), NUM_DEPENDENT_JOINTS);

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.col(i) = m.col(dependentJointIds[i]);
    }
}

void DigitDynamicsConstraints::get_independent_columns(MatX& r, const MatX& m) {
    assert(m.cols() == modelPtr_->nv);
    
    r.resize(m.rows(), NUM_INDEPENDENT_JOINTS);

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.col(i) = m.col(independentJointIds[i]);
    }
}

void DigitDynamicsConstraints::get_dependent_rows(MatX& r, const MatX& m) {
    assert(m.rows() == modelPtr_->nv);

    r.resize(NUM_DEPENDENT_JOINTS, m.cols());

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        r.row(i) = m.row(dependentJointIds[i]);
    }
}

void DigitDynamicsConstraints::get_independent_rows(MatX& r, const MatX& m) {
    assert(m.rows() == modelPtr_->nv);

    r.resize(NUM_INDEPENDENT_JOINTS, m.cols());

    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        r.row(i) = m.row(independentJointIds[i]);
    }
}

void DigitDynamicsConstraints::setupJointPosition(VecX& q, bool compute_derivatives) {
    if (recoverSavedData(q, compute_derivatives)) {
        return;
    }
    
    qcopy = q;

    // fill in dependent joint positions 
    // we first use approximation to give an initial guess for dependent joints
    // then we use gsl multidimensional root-finding to provide a more accurate solution

    // a Fourier based approximation. the code is modified from code generation by Matlab symbolic toolbox
    // right toe close loop
    double q1 = q(modelPtr_->getJointId("right_toe_A") - 1);
    double q2 = q(modelPtr_->getJointId("right_toe_B") - 1);
    double t2 = cos(q1);
    double t3 = cos(q2);
    double t4 = sin(q1);
    double t5 = sin(q2);
    double t6 = q1*2.0;
    double t7 = q2*2.0;
    double t8 = cos(t6);
    double t9 = cos(t7);
    double t10 = sin(t6);
    double t11 = sin(t7);
    qcopy(modelPtr_->getJointId("right_toe_A_rod") - 1) = -HigherOrderDerivatives::safeasin(t2*5.114967905744713E-1+t3*4.072902857734452E-1+t4*5.15019991000208E-1-t5*7.189052210514619E-1-t8*1.914308797309167E-1+t9*1.020482430826409E-1+t10*2.762037569148706E-1+t11*4.461078267300372E-1-t2*t3*5.669643819573642E-1+t2*t5*9.826313004441752E-1+t3*t4*7.066611122892024E-1-t4*t5*8.751268427354268E-2-t2*t9*1.526013089230614E-1+t3*t8*2.146475284950957E-1-t2*t11*6.373722738059316E-1-t3*t10*4.061233553924697E-1-t4*t9*2.591685269591633E-1-t5*t8*3.141536266927853E-1+t4*t11*1.154141052196728E-1+t5*t10*5.111686735370716E-2+t8*t9*2.629326893467769E-2+t8*t11*2.261999172185061E-1+t9*t10*1.714455893661597E-1-t10*t11*5.322298539488757E-2-3.496118766033006E-1);
    qcopy(modelPtr_->getJointId("right_A2") - 1) = HigherOrderDerivatives::safeasin(t2*8.08508988653757+t3*8.033204820099098+t4*1.820663506981432+t5*9.760946687487468E-1-t8*2.143120822774579-t9*1.998293355526279-t10*1.095961671709748-t11*4.953901287885355E-1-t2*t3*1.273036938982479E+1-t2*t5*1.117825053271672-t3*t4*2.420405429552186+t4*t5*1.170037217867019+t2*t9*2.835751824751641+t3*t8*3.021210259875125+t2*t11*5.805706025081291E-1+t3*t10*1.468854423683356+t4*t9*5.420356015832145E-1+t5*t8*2.326217949887465E-1+t4*t11*1.600391744429722E-1+t5*t10*1.170908991628156E-1-t8*t9*4.40337591097512E-1-t8*t11*1.065731420572547E-1-t9*t10*3.400548068856727E-1-t10*t11*4.121086576936894E-1-4.663948776295725);
    qcopy(modelPtr_->getJointId("right_toe_B_rod") - 1) = -HigherOrderDerivatives::safeasin(t2*(-6.277469604269263E-1)-t3*6.676022839778379E-1-t4*8.399423316628958E-1+t5*1.529719504893023E-1+t8*1.44381731007234E-1+t9*4.870768316848186E-1+t10*5.196004288358185E-1+t11*4.899787657102767E-1+t2*t3*1.038356653790214+t2*t5*1.224128324227442+t3*t4*1.136249345040226-t4*t5*2.962244825741616E-1-t2*t9*7.061936994304631E-1-t3*t8*2.508516890107657E-1-t2*t11*7.147313228957868E-1-t3*t10*7.361945978788478E-1-t4*t9*3.638713473685915E-1-t5*t8*4.137119943797418E-1+t4*t11*2.198247731883075E-1+t5*t10*1.190720324216782E-1+t8*t9*1.91374695742122E-1+t8*t11*2.674233531253601E-1+t9*t10*2.592073581826719E-1-t10*t11*1.033333326963321E-1+3.90361533934657E-1);
    qcopy(modelPtr_->getJointId("right_B2") - 1) = HigherOrderDerivatives::safeasin(t2*9.834460864042423+t3*1.012210878062529E+1+t4*4.114199768067971+t5*2.918615505121678-t8*2.33635398066763-t9*2.824506662141518-t10*2.436915399445716-t11*1.657161229950764-t2*t3*1.56824564483544E+1-t2*t5*4.071020737587724-t3*t4*5.741194151827515+t4*t5*1.412990674523749+t2*t9*3.933187185458805+t3*t8*3.362495208957303+t2*t11*2.232327688060054+t3*t10*3.377470435585635+t4*t9*1.500870022094839+t5*t8*1.104866050029796+t4*t11*1.111537011807542E-1+t5*t10*1.885808892070127E-1-t8*t9*5.75161230009884E-1-t8*t11*6.064451339600673E-1-t9*t10*9.059551372991537E-1-t10*t11*4.795029892364048E-1-5.829415122169579);
    qcopy(modelPtr_->getJointId("right_toe_pitch") - 1) = -HigherOrderDerivatives::safeasin(t2*1.893989559445235E+1+t3*1.920869550741533E+1+t4*6.868951098405778+t5*3.51561215609706-t8*4.080342460091422-t9*4.456343821978674-t10*3.407846365996183-t11*2.503419572692307-t2*t3*2.929974638198693E+1-t2*t5*5.157359318839895-t3*t4*8.730366732237334+t4*t5*1.967066798625298+t2*t9*5.895002722892592+t3*t8*5.481616764469075+t2*t11*3.209642327668312+t3*t10*4.724342470672166+t4*t9*2.235861962074205+t5*t8*1.242934577497499+t4*t11*1.119458771782081+t5*t10*1.139096075125015-t8*t9*3.306304563913454E-1-t8*t11*7.650487856200098E-1-t9*t10*1.240832663118515-t10*t11*1.549291982370939-1.135688467845833E+1);
    qcopy(modelPtr_->getJointId("right_toe_roll") - 1) = -HigherOrderDerivatives::safeasin(t2*1.925184467728597+t3*7.345366777287639+t4*2.843699508105022E+1+t5*2.720900626308319E+1+t8*3.57433462778918-t9*6.332457856504018-t10*1.691966583575934E+1-t11*1.60681725794845E+1-t2*t3*6.299462363896048-t2*t5*4.117199083152367E+1-t3*t4*4.282360516873334E+1+t4*t5*7.910198316638063E-2+t2*t9*7.854243405249483-t3*t8*4.011841131466795+t2*t11*2.320584908063003E+1+t3*t10*2.435545111239733E+1+t4*t9*1.284334673751474E+1+t5*t8*1.242519321726289E+1-t4*t11*1.106503053903957+t5*t10*8.20592115289103E-1-t8*t9*6.089534825702931E-1-t8*t11*7.206363286390523-t9*t10*7.502224690460121+t10*t11*1.255468934813467E-1-3.446675672448973);

    // left toe close loop
    q1 = q(modelPtr_->getJointId("left_toe_A") - 1);
    q2 = q(modelPtr_->getJointId("left_toe_B") - 1);
    t2 = cos(q1);
    t3 = cos(q2);
    t4 = sin(q1);
    t5 = sin(q2);
    t6 = q1*2.0;
    t7 = q2*2.0;
    t8 = cos(t6);
    t9 = cos(t7);
    t10 = sin(t6);
    t11 = sin(t7);
    qcopy(modelPtr_->getJointId("left_toe_A_rod") - 1) = -HigherOrderDerivatives::safeasin(t2*(-5.114926838701811E-1)-t3*4.072861778028993E-1+t4*5.150198314974532E-1-t5*7.189053388109297E-1+t8*1.914297694169048E-1-t9*1.020493568564005E-1+t10*2.762038362432678E-1+t11*4.461078824742524E-1+t2*t3*5.669587551620541E-1+t2*t5*9.826314482751646E-1+t3*t4*7.066613185154341E-1+t4*t5*8.751254943270917E-2+t2*t9*1.526028367712922E-1-t3*t8*2.146460048968412E-1-t2*t11*6.3737234338624E-1-t3*t10*4.061234579946822E-1-t4*t9*2.591685735898548E-1-t5*t8*3.141536565807085E-1-t4*t11*1.15414032157577E-1-t5*t10*5.111679397817798E-2-t8*t9*2.629368451925844E-2+t8*t11*2.261999309581469E-1+t9*t10*1.714456126031132E-1+t10*t11*5.322294548977349E-2+3.49608876912165E-1);
    qcopy(modelPtr_->getJointId("left_A2") - 1) = HigherOrderDerivatives::safeasin(t2*8.085018160207163+t3*8.033133064086231-t4*1.820662389025083-t5*9.760941431166779E-1-t8*2.143101587807895-t9*1.998274051537184+t10*1.095961134647075+t11*4.953899201857669E-1-t2*t3*1.273027136495291E+1+t2*t5*1.11782444127773+t3*t4*2.420403987960523+t4*t5*1.17003991016308+t2*t9*2.835725404415498+t3*t8*3.021183922843985-t2*t11*5.805703711334091E-1-t3*t10*1.468853731639082-t4*t9*5.420352775054091E-1-t5*t8*2.326217084069311E-1+t4*t11*1.600377315075955E-1+t5*t10*1.170894522681423E-1-t8*t9*4.403304524513951E-1+t8*t11*1.065731191784157E-1+t9*t10*3.400546515442943E-1-t10*t11*4.121078791939159E-1-4.663896238435751);
    qcopy(modelPtr_->getJointId("left_toe_B_rod") - 1) = -HigherOrderDerivatives::safeasin(t2*6.277394162562625E-1+t3*6.6759473381467E-1-t4*8.399421448959555E-1+t5*1.52972128181329E-1-t8*1.443797201746081E-1-t9*4.870748185116889E-1+t10*5.196003431996498E-1+t11*4.899786813403554E-1-t2*t3*1.038346357873659+t2*t5*1.224128102824049+t3*t4*1.13624910964581+t4*t5*2.962247949009167E-1+t2*t9*7.06190949266457E-1+t3*t8*2.508489399062072E-1-t2*t11*7.147312184882999E-1-t3*t10*7.361944907642444E-1-t4*t9*3.638712999144287E-1-t5*t8*4.137119506580644E-1-t4*t11*2.198249429188003E-1-t5*t10*1.190722007592886E-1-t8*t9*1.913739574570494E-1+t8*t11*2.67423333074767E-1+t9*t10*2.592073373292207E-1+t10*t11*1.033334244581574E-1-3.903559984881447E-1);
    qcopy(modelPtr_->getJointId("left_B2") - 1) = HigherOrderDerivatives::safeasin(t2*9.83436298937519+t3*1.012201086685809E+1-t4*4.114198671144258-t5*2.918615149671792-t8*2.3363277301986-t9*2.824480323210008+t10*2.436914889769644+t11*1.657161129062592-t2*t3*1.568232268107518E+1+t2*t5*4.071020370775729+t3*t4*5.741192746046224+t4*t5*1.412994342972937+t2*t9*3.933151135492005+t3*t8*3.362459265404445-t2*t11*2.232327610328257-t3*t10*3.377469784128673-t4*t9*1.500869712732554-t5*t8*1.104866037480153+t4*t11*1.111517321363341E-1+t5*t10*1.885789156342841E-1-t8*t9*5.751514893646594E-1+t8*t11*6.064451565041576E-1+t9*t10*9.05954995114355E-1-t10*t11*4.795019258370027E-1-5.829343436697216);
    qcopy(modelPtr_->getJointId("left_toe_pitch") - 1) = HigherOrderDerivatives::safeasin(t2*1.893971200186481E+1+t3*1.920851184673561E+1-t4*6.868948549897667-t5*3.515611050123151-t8*4.080293186785779-t9*4.456294384490135+t10*3.407845158508307+t11*2.503419165594914-t2*t3*2.92994954151725E+1+t2*t5*5.157358065307949+t3*t4*8.730363458834317+t4*t5*1.967073609987012+t2*t9*5.894935044834759+t3*t8*5.481549285471837-t2*t11*3.209641901519396-t3*t10*4.72434092223045-t4*t9*2.235861236273534-t5*t8*1.242934428625482+t4*t11*1.119455114838831+t5*t10*1.139092408157223-t8*t9*3.306121591996909E-1+t8*t11*7.650487658943779E-1+t9*t10*1.240832321393588-t10*t11*1.549290006003741-1.135675024134487E+1);
    qcopy(modelPtr_->getJointId("left_toe_roll") - 1) = HigherOrderDerivatives::safeasin(t2*1.92500223858939+t3*7.345184505784512-t4*2.843699897601034E+1-t5*2.720901076343499E+1+t8*3.574383636815427-t9*6.332408732694388+t10*1.691966794969343E+1+t11*1.606817500398276E+1-t2*t3*6.299213096460508+t2*t5*4.117199681273369E+1+t3*t4*4.282361029947926E+1+t4*t5*7.910853681196074E-2+t2*t9*7.854176116345911-t3*t8*4.011908279005063-t2*t11*2.320585230931768E+1-t3*t10*2.435545390288711E+1-t4*t9*1.284334796900915E+1-t5*t8*1.242519468406818E+1-t4*t11*1.106506602910529+t5*t10*8.205885613072367E-1-t8*t9*6.089352667535674E-1+t8*t11*7.206364083227697+t9*t10*7.502225364895965+t10*t11*1.255488247882333E-1-3.446542349892159);

    // right knee close loop
    q1 = q(modelPtr_->getJointId("right_knee") - 1);
    t2 = cos(q1);
    t3 = sin(q1);
    t4 = q1*2.0;
    t5 = cos(t4);
    t6 = sin(t4);
    qcopy(modelPtr_->getJointId("right_tarsus") - 1) = -HigherOrderDerivatives::safeasin(t2*1.155848969647063E-3+t3*1.004686948291003+t5*1.274417498011625E-4-t6*1.785981355062532E-3-1.132590494159057E-2);
    qcopy(modelPtr_->getJointId("right_achilles_rod") - 1) = -HigherOrderDerivatives::safeasin(t2*(-1.587289102030986E-3)-t3*1.001736520672665+t5*3.407131509821247E-4+t6*9.646678263881318E-4+1.539911054998293E-3);
    qcopy(modelPtr_->getJointId("right_ach2") - 1) = HigherOrderDerivatives::safeasin(t2*(-7.197863326636346E-2)-t3*8.929579539511397E-3+t5*2.669904889172627E-4+t6*8.46571305589265E-5+7.18964949007849E-2);

    // left knee close loop
    q1 = q(modelPtr_->getJointId("left_knee") - 1);
    t2 = cos(q1);
    t3 = sin(q1);
    t4 = q1*2.0;
    t5 = cos(t4);
    t6 = sin(t4);
    qcopy(modelPtr_->getJointId("left_tarsus") - 1) = HigherOrderDerivatives::safeasin(t2*1.155848972188414E-3-t3*1.004686948291033+t5*1.274417489907877E-4+t6*1.785981355072367E-3-1.132590494335349E-2);
    qcopy(modelPtr_->getJointId("left_achilles_rod") - 1) = HigherOrderDerivatives::safeasin(t2*(-1.587289102219775E-3)+t3*1.001736520672708+t5*3.407131510426615E-4-t6*9.646678264174077E-4+1.539911055129276E-3);
    qcopy(modelPtr_->getJointId("left_ach2") - 1) = HigherOrderDerivatives::safeasin(t2*(-7.197863326651638E-2)+t3*8.929579539517018E-3+t5*2.669904889661342E-4-t6*8.465713056222183E-5+7.189649490089099E-2);

    fkPtr_->compute(modelPtr_->getJointId("Rz"), 
                    contact_joint_id, 
                    qcopy, 
                    nullptr, 
                    &stance_foot_endT, 
                    0);
    Transform torso_T = stance_foot_T_des * fkPtr_->getTransform().inverse();
    qcopy.head(6) = torso_T.getXYZRPY();

    // gsl multidimensional root-finding
    const gsl_multiroot_fdfsolver_type *T;
    gsl_multiroot_fdfsolver *s;

    int status;
    size_t i, iter = 0;

    const size_t n = NUM_DEPENDENT_JOINTS;
    gsl_multiroot_function_fdf f = {&fillDependent_f,
                                    &fillDependent_df,
                                    &fillDependent_fdf,
                                    n, this};

    gsl_vector *x = gsl_vector_alloc(n);

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        gsl_vector_set(x, i, qcopy(dependentJointIds[i]));
    }

    T = gsl_multiroot_fdfsolver_hybridsj;
    s = gsl_multiroot_fdfsolver_alloc(T, n);
    gsl_multiroot_fdfsolver_set(s, &f, x);

    // printf("\nstart gsl iterations:\n");

    do {
        iter++;
        status = gsl_multiroot_fdfsolver_iterate(s);

        if (status) break;

        status = gsl_multiroot_test_residual(s->f, 1e-14);

        // printf ("iter = %ld, status = %s\n", iter, gsl_strerror(status));
        // for (int i = 0; i < 6; i++) {
        //     printf("%f, ", gsl_vector_get(s->x, i)) ;
        // }
        // for (int i = 0; i < 6; i++) {
        //     printf("%f, ", gsl_vector_get(s->f, 18 + i)) ;
        // }
        // printf("\n");
    }
    while (status == GSL_CONTINUE && iter < 50);

    // printf ("total iter = %ld, status = %s\n\n", iter, gsl_strerror(status));

    // the optimal solution found by gsl!
    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        qcopy(dependentJointIds[i]) = gsl_vector_get(s->x, i);
    }

    gsl_multiroot_fdfsolver_free(s);
    gsl_vector_free(x);

    q = qcopy;

    if (compute_derivatives) {
        get_c(q);
        get_J(q);

        get_dependent_columns(J_dep, J);
        get_independent_columns(J_indep, J);

        J_dep_qr = QRSolver(J_dep);
        J_dep_T_qr = QRSolver(J_dep.transpose());

        // sanity check on uniqueness (these two arguments are actually equivalent)
        assert(J_dep_qr.rank() == J_dep.rows() && J_dep_qr.rank() == J_dep.cols());
        assert(J_dep_T_qr.rank() == J_dep.rows() && J_dep_T_qr.rank() == J_dep.cols());

        P_dep = -J_dep_qr.solve(J_indep);
        pq_dep_pq_indep = P_dep;
    }
    
    updateQueue(q, compute_derivatives);
}

int fillDependent_f(const gsl_vector* x, void *params, gsl_vector* f) {
    DigitDynamicsConstraints* constraintsData = (DigitDynamicsConstraints*)params;

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        constraintsData->qcopy(constraintsData->dependentJointIds[i]) = gsl_vector_get(x, i);
    }

    constraintsData->get_c(constraintsData->qcopy);

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        gsl_vector_set(f, i, constraintsData->c(i));
    }

    return GSL_SUCCESS;
}

int fillDependent_df(const gsl_vector* x, void *params, gsl_matrix* J) {
    DigitDynamicsConstraints* constraintsData = (DigitDynamicsConstraints*)params;

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        constraintsData->qcopy(constraintsData->dependentJointIds[i]) = gsl_vector_get(x, i);
    }

    constraintsData->get_J(constraintsData->qcopy);

    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        for (int j = 0; j < NUM_DEPENDENT_JOINTS; j++) {
            gsl_matrix_set(J, i, j, constraintsData->J(i, constraintsData->dependentJointIds[j]));
        }
    }

    return GSL_SUCCESS;
}

int fillDependent_fdf(const gsl_vector* x, void *params, gsl_vector* f, gsl_matrix* J) {
    fillDependent_f(x, params, f);
    fillDependent_df(x, params, J);

    return GSL_SUCCESS;
}

void DigitDynamicsConstraints::get_c(const VecX& q) {
    double t2 = cos(q(9));
    double t3 = cos(q(10));
    double t4 = cos(q(11));
    double t5 = cos(q(12));
    double t6 = cos(q(13));
    double t7 = cos(q(14));
    double t8 = cos(q(15));
    double t9 = cos(q(16));
    double t10 = cos(q(17));
    double t11 = cos(q(18));
    double t12 = cos(q(19));
    double t13 = cos(q(20));
    double t14 = cos(q(24));
    double t15 = cos(q(25));
    double t16 = cos(q(26));
    double t17 = cos(q(27));
    double t18 = cos(q(28));
    double t19 = cos(q(29));
    double t20 = cos(q(30));
    double t21 = cos(q(31));
    double t22 = cos(q(32));
    double t23 = cos(q(33));
    double t24 = cos(q(34));
    double t25 = cos(q(35));
    double t26 = sin(q(9));
    double t27 = sin(q(10));
    double t28 = sin(q(11));
    double t29 = sin(q(12));
    double t30 = sin(q(13));
    double t31 = sin(q(14));
    double t32 = sin(q(15));
    double t33 = sin(q(16));
    double t34 = sin(q(17));
    double t35 = sin(q(18));
    double t36 = sin(q(19));
    double t37 = sin(q(20));
    double t38 = sin(q(24));
    double t39 = sin(q(25));
    double t40 = sin(q(26));
    double t41 = sin(q(27));
    double t42 = sin(q(28));
    double t43 = sin(q(29));
    double t44 = sin(q(30));
    double t45 = sin(q(31));
    double t46 = sin(q(32));
    double t47 = sin(q(33));
    double t48 = sin(q(34));
    double t49 = sin(q(35));
    double t50 = t2*t2;
    double t51 = t3*t3;
    double t52 = t4*t4;
    double t53 = t5*t5;
    double t54 = t6*t6;
    double t55 = t7*t7;
    double t56 = t8*t8;
    double t57 = t9*t9;
    double t58 = t10*t10;
    double t59 = t11*t11;
    double t60 = t12*t12;
    double t61 = t13*t13;
    double t62 = t14*t14;
    double t63 = t15*t15;
    double t64 = t16*t16;
    double t65 = t17*t17;
    double t66 = t26*t26;
    double t67 = t27*t27;
    double t68 = t28*t28;
    double t69 = t29*t29;
    double t70 = t30*t30;
    double t71 = t31*t31;
    double t72 = t32*t32;
    double t73 = t33*t33;
    double t74 = t34*t34;
    double t75 = t35*t35;
    double t76 = t36*t36;
    double t77 = t37*t37;
    double t78 = t38*t38;
    double t79 = t39*t39;
    double t80 = t40*t40;
    double t81 = t41*t41;
    double t82 = t24*1.696216709330505E-2;
    double t83 = t48*1.696216709330505E-2;
    double t84 = t49*-9.551E-3;
    double t85 = t49*9.551E-3;
    double t86 = t24*-5.143951577823025E-2;
    double t87 = t24*5.143951577823025E-2;
    double t88 = t48*5.143951577823025E-2;
    double t97 = t12*t36*-1.758834080497775E-1;
    double t98 = t12*t36*1.758834080497775E-1;
    double t99 = t12*t36*-1.052664560664766E-1;
    double t100 = t12*t36*1.052664560664766E-1;
    double t89 = t76*-4.263322803323831E-2;
    double t90 = t76*4.263322803323831E-2;
    double t91 = t61*-2.04E-1;
    double t92 = t61*2.04E-1;
    double t93 = t77*-2.04E-1;
    double t94 = t77*2.04E-1;
    double t95 = t76*-1.405829597511123E-2;
    double t96 = t76*1.405829597511123E-2;
    double t101 = t60*6.263322803323831E-2;
    double t102 = t61*2.0E-2;
    double t103 = t60*-1.899417040248888E-1;
    double t104 = t60*1.899417040248888E-1;
    double t105 = t77*2.0E-2;
    double t106 = t12*t61*1.696216709330505E-2;
    double t107 = t12*t77*1.696216709330505E-2;
    double t108 = t36*t61*1.696216709330505E-2;
    double t109 = t37*t60*9.551E-3;
    double t110 = t36*t77*1.696216709330505E-2;
    double t111 = t12*t61*5.143951577823025E-2;
    double t112 = t37*t76*9.551E-3;
    double t113 = t12*t77*5.143951577823025E-2;
    double t114 = t36*t61*-5.143951577823025E-2;
    double t115 = t36*t61*5.143951577823025E-2;
    double t116 = t36*t77*-5.143951577823025E-2;
    double t117 = t36*t77*5.143951577823025E-2;
    double t118 = t61*t76*6.263322803323831E-2;
    double t119 = t60*t61*-4.263322803323831E-2;
    double t120 = t60*t61*4.263322803323831E-2;
    double t121 = t76*t77*6.263322803323831E-2;
    double t122 = t60*t77*-4.263322803323831E-2;
    double t123 = t60*t77*4.263322803323831E-2;
    double t124 = t61*t76*-1.899417040248888E-1;
    double t125 = t61*t76*1.899417040248888E-1;
    double t126 = t60*t61*-1.405829597511123E-2;
    double t127 = t60*t61*1.405829597511123E-2;
    double t128 = t76*t77*-1.899417040248888E-1;
    double t129 = t76*t77*1.899417040248888E-1;
    double t130 = t60*t77*-1.405829597511123E-2;
    double t131 = t60*t77*1.405829597511123E-2;
    double t132 = t61*t98;
    double t133 = t77*t98;
    double t134 = t61*t100;
    double t135 = t77*t100;

    c[0] = t18*1.034589188110661E-3+t42*5.699060997402858E-2+t82+t88-t24*t25*9.070578524442013E-3-t18*t44*2.246221860400801E-3-t24*t49*1.69996184260823E-2+t25*t48*2.991020934719675E-3-t42*t44*1.26575449899492E-2+t48*t49*5.605619802270151E-3+t18*t19*t20*3.399783924207052E-1+t18*t20*t43*3.105990081579729E-3+t19*t20*t42*3.020283789547073E-3-t20*t42*t43*3.397508858570615E-1-3.49E-1;
    c[1] = t18*5.699060997402858E-2-t42*1.034589188110661E-3+t83+t86-t24*t25*2.991020934719675E-3-t18*t44*1.26575449899492E-2-t24*t49*5.605619802270151E-3-t25*t48*9.070578524442013E-3+t42*t44*2.246221860400801E-3-t48*t49*1.69996184260823E-2+t18*t19*t20*3.020283789547073E-3-t18*t20*t43*3.397508858570615E-1-t19*t20*t42*3.399783924207052E-1-t20*t42*t43*3.105990081579729E-3-6.0E-3;
    c[2] = t25*1.79E-2+t44*3.39756885202024E-1+t84+t19*t20*2.360206106172353E-3-t20*t43*1.263678697118524E-2-1.96E-2;
    c[3] = t21*1.034589188110661E-3-t45*5.699060997402856E-2+t82+t88-t24*t25*9.070578524442012E-3-t21*t47*8.234792310892687E-4+t24*t49*1.718955829676478E-2+t25*t48*2.991020934719677E-3-t45*t47*1.549180159908108E-2-t48*t49*5.668252425759201E-3+t21*t22*t23*2.879825211612492E-1-t21*t23*t46*3.064210757541298E-3-t22*t23*t45*3.104080344633556E-3-t23*t45*t46*2.87566285868891E-1-2.97E-1;
    c[4] = t21*5.699060997402856E-2+t45*1.034589188110661E-3+t83+t86-t24*t25*2.991020934719677E-3+t21*t47*1.549180159908108E-2+t24*t49*5.668252425759201E-3-t25*t48*9.070578524442012E-3-t45*t47*8.234792310892687E-4+t48*t49*1.718955829676478E-2+t21*t22*t23*3.104080344633556E-3+t21*t23*t46*2.87566285868891E-1+t22*t23*t45*2.879825211612492E-1-t23*t45*t46*3.064210757541298E-3-6.0E-3;
    c[5] = t25*-1.81E-2-t47*2.875818595898751E-1+t84-t22*t23*6.574122182670539E-4+t23*t46*1.549973690114125E-2+1.96E-2;
    c[6] = t55*4.075092868664323E-5+t56*2.945782751508952E-2+t71*1.421556223845944E-6+t72*2.945782751508952E-2+t91+t93+t95+t99+t103+t106+t107+t114+t116+t124+t126+t128+t130+t134+t135-t12*t13*9.070578524442013E-3-t7*t31*1.522231734027378E-5+t12*t37*1.69996184260823E-2-t13*t36*2.991020934719675E-3+t6*t55*7.74045457557333E-4+t6*t56*4.455658896021977E-4+t36*t37*5.605619802270151E-3-t6*t71*1.8502215904887E-4+t6*t72*4.455658896021977E-4-t30*t55*4.104331657223103E-4-t30*t56*2.809111102928818E-2-t30*t71*2.848906577901809E-2-t30*t72*2.809111102928818E-2+t54*t55*2.961197292377137E-2+t54*t56*3.670379334125467E-5+t55*t56*1.421556223845944E-6-t54*t71*1.486767171126281E-4-t55*t70*1.527238524580169E-4+t54*t72*3.670379334125467E-5+t56*t70*5.468691569234733E-6+t55*t72*1.421556223845944E-6+t56*t71*4.075092868664323E-5+t70*t71*2.964725516088879E-2+t70*t72*5.468691569234733E-6+t71*t72*4.075092868664323E-5+t6*t7*t8*3.399783924207052E-1+t6*t7*t31*2.809726196867663E-2-t6*t8*t31*3.105990081579729E-3-t7*t8*t30*3.020283789547073E-3-t7*t30*t31*9.515128902424641E-4-t8*t30*t31*3.397508858570615E-1+t6*t30*t55*1.643509328293678E-2+t6*t30*t56*3.732977901401361E-5+t7*t31*t54*1.646136108951249E-2-t6*t32*t55*2.246221860400801E-3+t7*t31*t56*1.522231734027378E-5-t6*t30*t71*1.64724230619508E-2+t6*t30*t72*3.732977901401361E-5-t7*t31*t70*1.644613877217222E-2-t6*t32*t71*2.246221860400801E-3+t7*t31*t72*1.522231734027378E-5-t6*t55*t56*1.8502215904887E-4+t30*t32*t55*1.26575449899492E-2-t6*t55*t72*1.8502215904887E-4+t6*t56*t71*7.74045457557333E-4+t30*t32*t71*1.26575449899492E-2-t30*t55*t56*2.848906577901809E-2+t6*t71*t72*7.74045457557333E-4-t30*t55*t72*2.848906577901809E-2-t30*t56*t71*4.104331657223103E-4-t54*t55*t56*1.486767171126281E-4-t30*t71*t72*4.104331657223103E-4-t54*t55*t72*1.486767171126281E-4+t54*t56*t71*2.961197292377137E-2+t55*t56*t70*2.964725516088879E-2+t54*t71*t72*2.961197292377137E-2+t55*t70*t72*2.964725516088879E-2-t56*t70*t71*1.527238524580169E-4-t70*t71*t72*1.527238524580169E-4-t6*t7*t30*t31*5.956062305521773E-2-t6*t7*t31*t56*2.809726196867663E-2-t6*t7*t31*t72*2.809726196867663E-2+t7*t30*t31*t56*9.515128902424641E-4+t7*t30*t31*t72*9.515128902424641E-4-t6*t30*t55*t56*1.64724230619508E-2-t7*t31*t54*t56*1.646136108951249E-2-t6*t30*t55*t72*1.64724230619508E-2+t6*t30*t56*t71*1.643509328293678E-2-t7*t31*t54*t72*1.646136108951249E-2+t7*t31*t56*t70*1.644613877217222E-2+t6*t30*t71*t72*1.643509328293678E-2+t7*t31*t70*t72*1.644613877217222E-2+t6*t7*t30*t31*t56*5.956062305521773E-2+t6*t7*t30*t31*t72*5.956062305521773E-2;
    c[7] = t55*-2.348358602281135E-5-t56*1.697569721208548E-2-t71*8.192018917078324E-7-t72*1.697569721208548E-2+t89+t97+t101+t102+t105+t108+t110+t111+t113+t118+t119+t121+t122+t132+t133+t12*t13*2.991020934719675E-3+t7*t31*8.772182874056074E-6-t12*t37*5.605619802270151E-3-t13*t36*9.070578524442013E-3-t6*t55*4.104331657223105E-4-t6*t56*2.809111102928818E-2+t36*t37*1.69996184260823E-2-t6*t71*2.848906577901809E-2-t6*t72*2.809111102928818E-2-t30*t55*7.74045457557333E-4-t30*t56*4.455658896021977E-4+t30*t71*1.8502215904887E-4-t30*t72*4.455658896021977E-4-t54*t55*2.707115655202057E-4+t54*t56*6.513495549747348E-6-t55*t56*8.192018917078324E-7-t54*t71*1.673580193002955E-2-t55*t70*1.670580484845698E-2+t54*t72*6.513495549747348E-6-t56*t70*3.081628346426626E-5-t55*t72*8.192018917078324E-7-t56*t71*2.348358602281135E-5-t70*t71*2.633788680787505E-4-t70*t72*3.081628346426626E-5-t71*t72*2.348358602281135E-5-t6*t7*t8*3.020283789547073E-3-t6*t7*t31*9.51512890242464E-4-t6*t8*t31*3.397508858570615E-1-t7*t8*t30*3.399783924207052E-1-t7*t30*t31*2.809726196867663E-2+t8*t30*t31*3.105990081579729E-3-t6*t30*t55*2.976469677622939E-2-t6*t30*t56*3.123510177201994E-5-t7*t31*t54*2.978469761904589E-2+t6*t32*t55*1.26575449899492E-2-t7*t31*t56*8.772182874056074E-6+t6*t30*t71*2.979593187800142E-2-t6*t30*t72*3.123510177201994E-5+t7*t31*t70*2.977592543617184E-2+t6*t32*t71*1.26575449899492E-2-t7*t31*t72*8.772182874056074E-6-t6*t55*t56*2.848906577901809E-2+t30*t32*t55*2.246221860400801E-3-t6*t55*t72*2.848906577901809E-2-t6*t56*t71*4.104331657223105E-4+t30*t32*t71*2.246221860400801E-3+t30*t55*t56*1.8502215904887E-4-t6*t71*t72*4.104331657223105E-4+t30*t55*t72*1.8502215904887E-4-t30*t56*t71*7.74045457557333E-4-t54*t55*t56*1.673580193002955E-2-t30*t71*t72*7.74045457557333E-4-t54*t55*t72*1.673580193002955E-2-t54*t56*t71*2.707115655202057E-4-t55*t56*t70*2.633788680787505E-4-t54*t71*t72*2.707115655202057E-4-t55*t70*t72*2.633788680787505E-4-t56*t70*t71*1.670580484845698E-2-t70*t71*t72*1.670580484845698E-2-t6*t7*t30*t31*3.290749986168471E-2+t6*t7*t31*t56*9.51512890242464E-4+t6*t7*t31*t72*9.51512890242464E-4+t7*t30*t31*t56*2.809726196867663E-2+t7*t30*t31*t72*2.809726196867663E-2+t6*t30*t55*t56*2.979593187800142E-2+t7*t31*t54*t56*2.978469761904589E-2+t6*t30*t55*t72*2.979593187800142E-2-t6*t30*t56*t71*2.976469677622939E-2+t7*t31*t54*t72*2.978469761904589E-2-t7*t31*t56*t70*2.977592543617184E-2-t6*t30*t71*t72*2.976469677622939E-2-t7*t31*t70*t72*2.977592543617184E-2+t6*t7*t30*t31*t56*3.290749986168471E-2+t6*t7*t30*t31*t72*3.290749986168471E-2;
    c[8] = t109+t112-t6*t55*6.213602549397351E-4+t6*t56*8.271781369638137E-4+t13*t60*1.79E-2-t6*t71*2.058178820240719E-4+t6*t72*8.271781369638137E-4+t30*t55*1.101395785003593E-3-t30*t56*9.852121016062406E-4+t13*t76*1.79E-2-t30*t71*1.161836833973452E-4-t30*t72*9.852121016062406E-4-t54*t55*1.084459637073036E-2+t54*t56*1.04947201084778E-3-t54*t71*9.804875640117422E-3-t55*t70*1.084459637073035E-2+t54*t72*1.04947201084778E-3+t56*t70*1.04947201084778E-3-t70*t71*9.804875640117427E-3+t70*t72*1.04947201084778E-3+t6*t7*t31*1.218023269812134E-3+t7*t30*t31*4.163488242625716E-4+t7*t8*t54*2.360206106172353E-3+t7*t8*t70*2.360206106172353E-3+t7*t31*t54*2.212067055639549E-4+t8*t31*t54*1.263678697118524E-2+t7*t31*t70*2.212067055639517E-4+t8*t31*t70*1.263678697118524E-2-t6*t55*t56*2.058178820240719E-4-t6*t55*t72*2.058178820240719E-4-t6*t56*t71*6.213602549397351E-4-t30*t55*t56*1.161836833973452E-4+t32*t54*t55*3.39756885202024E-1-t6*t71*t72*6.213602549397351E-4-t30*t55*t72*1.161836833973452E-4+t30*t56*t71*1.101395785003593E-3+t32*t54*t71*3.39756885202024E-1+t32*t55*t70*3.39756885202024E-1-t54*t55*t56*9.804875640117422E-3+t30*t71*t72*1.101395785003593E-3+t32*t70*t71*3.39756885202024E-1-t54*t55*t72*9.804875640117422E-3-t54*t56*t71*1.084459637073036E-2-t55*t56*t70*9.804875640117427E-3-t54*t71*t72*1.084459637073036E-2-t55*t70*t72*9.804875640117427E-3-t56*t70*t71*1.084459637073035E-2-t70*t71*t72*1.084459637073035E-2-t6*t7*t31*t56*1.218023269812134E-3-t6*t7*t31*t72*1.218023269812134E-3-t7*t30*t31*t56*4.163488242625716E-4-t7*t30*t31*t72*4.163488242625716E-4-t7*t31*t54*t56*2.212067055639549E-4-t7*t31*t54*t72*2.212067055639549E-4-t7*t31*t56*t70*2.212067055639517E-4-t7*t31*t70*t72*2.212067055639517E-4;
    c[9] = t58*1.607521019272677E-4+t59*5.533895870788692E-2+t74*2.891901858161877E-7+t75*5.533895870788692E-2+t91+t93+t95+t99+t103+t106+t107+t114+t116+t124+t126+t128+t130+t134+t135-t12*t13*9.070578524442012E-3-t10*t34*1.363641158467861E-5-t12*t37*1.718955829676478E-2-t13*t36*2.991020934719677E-3+t9*t58*8.255702643738717E-4+t9*t59*4.936925916256938E-4-t36*t37*5.668252425759201E-3-t9*t74*2.846736678889049E-4+t9*t75*4.936925916256938E-4-t33*t58*4.353713116283091E-4+t33*t59*2.893932048270239E-2+t33*t74*2.848666080295449E-2+t33*t75*2.893932048270239E-2+t57*t58*5.551356637805161E-2+t57*t59*1.632022254576177E-4+t58*t59*2.891901858161877E-7-t57*t74*1.767686035092276E-4-t58*t73*1.743184799788775E-4+t57*t75*1.632022254576177E-4-t59*t73*2.160933344533858E-6+t58*t75*2.891901858161877E-7+t59*t74*1.607521019272677E-4+t73*t74*5.567647941332341E-2-t73*t75*2.160933344533858E-6+t74*t75*1.607521019272677E-4+t9*t10*t11*2.879825211612492E-1-t9*t10*t34*2.896401859571867E-2+t9*t11*t34*3.064210757541298E-3+t10*t11*t33*3.104080344633556E-3-t10*t33*t34*1.113095897646181E-3-t11*t33*t34*2.87566285868891E-1-t9*t33*t58*1.576769239528801E-2-t9*t33*t59*3.197767103277491E-5-t10*t34*t57*1.577685653751059E-2-t9*t35*t58*8.234792310892687E-4+t10*t34*t59*1.363641158467861E-5+t9*t33*t74*1.579967006632079E-2-t9*t33*t75*3.197767103277491E-5+t10*t34*t73*1.579049294909527E-2-t9*t35*t74*8.234792310892687E-4+t10*t34*t75*1.363641158467861E-5-t9*t58*t59*2.846736678889049E-4+t33*t35*t58*1.549180159908108E-2-t9*t58*t75*2.846736678889049E-4+t9*t59*t74*8.255702643738717E-4+t33*t35*t74*1.549180159908108E-2+t33*t58*t59*2.848666080295449E-2+t9*t74*t75*8.255702643738717E-4+t33*t58*t75*2.848666080295449E-2-t33*t59*t74*4.353713116283091E-4-t57*t58*t59*1.767686035092276E-4-t33*t74*t75*4.353713116283091E-4-t57*t58*t75*1.767686035092276E-4+t57*t59*t74*5.551356637805161E-2+t58*t59*t73*5.567647941332341E-2+t57*t74*t75*5.551356637805161E-2+t58*t73*t75*5.567647941332341E-2-t59*t73*t74*1.743184799788775E-4-t73*t74*t75*1.743184799788775E-4-t9*t10*t33*t34*1.115410112085772E-1+t9*t10*t34*t59*2.896401859571867E-2+t9*t10*t34*t75*2.896401859571867E-2+t10*t33*t34*t59*1.113095897646181E-3+t10*t33*t34*t75*1.113095897646181E-3+t9*t33*t58*t59*1.579967006632079E-2+t10*t34*t57*t59*1.577685653751059E-2+t9*t33*t58*t75*1.579967006632079E-2-t9*t33*t59*t74*1.576769239528801E-2+t10*t34*t57*t75*1.577685653751059E-2-t10*t34*t59*t73*1.579049294909527E-2-t9*t33*t74*t75*1.576769239528801E-2-t10*t34*t73*t75*1.579049294909527E-2+t9*t10*t33*t34*t59*1.115410112085772E-1+t9*t10*t33*t34*t75*1.115410112085772E-1;
    c[10] = t58*-4.923938257231622E-5-t59*1.695067203665005E-2-t74*8.85807776373908E-8-t75*1.695067203665005E-2+t89+t97+t101+t102+t105+t108+t110+t111+t113+t118+t119+t121+t122+t132+t133+t12*t13*2.991020934719677E-3+t10*t34*4.176918863775431E-6+t12*t37*5.668252425759201E-3-t13*t36*9.070578524442012E-3+t9*t58*4.353713116283091E-4-t9*t59*2.893932048270239E-2-t36*t37*1.718955829676478E-2-t9*t74*2.848666080295449E-2-t9*t75*2.893932048270239E-2+t33*t58*8.255702643738717E-4+t33*t59*4.936925916256938E-4-t33*t74*2.846736678889049E-4+t33*t75*4.936925916256938E-4-t57*t58*5.915341110698354E-4-t57*t59*8.675146158589342E-6-t58*t59*8.85807776373908E-8-t57*t74*1.639979074277157E-2-t58*t73*1.635922650635785E-2-t57*t75*8.675146158589342E-6-t59*t73*4.065281719136426E-5-t58*t75*8.85807776373908E-8-t59*t74*4.923938257231622E-5-t73*t74*6.001206764507874E-4-t73*t75*4.065281719136426E-5-t74*t75*4.923938257231622E-5-t9*t10*t11*3.104080344633556E-3+t9*t10*t34*1.113095897646181E-3+t9*t11*t34*2.87566285868891E-1+t10*t11*t33*2.879825211612492E-1-t10*t33*t34*2.896401859571867E-2+t11*t33*t34*3.064210757541298E-3+t9*t33*t58*5.568788485803049E-2+t9*t33*t59*1.653631588021515E-4+t10*t34*t57*5.576841714485674E-2-t9*t35*t58*1.549180159908108E-2-t10*t34*t59*4.176918863775431E-6-t9*t33*t74*5.585324801683264E-2+t9*t33*t75*1.653631588021515E-4-t10*t34*t73*5.577259406372051E-2-t9*t35*t74*1.549180159908108E-2-t10*t34*t75*4.176918863775431E-6-t9*t58*t59*2.848666080295449E-2-t33*t35*t58*8.234792310892687E-4-t9*t58*t75*2.848666080295449E-2+t9*t59*t74*4.353713116283091E-4-t33*t35*t74*8.234792310892687E-4-t33*t58*t59*2.846736678889049E-4+t9*t74*t75*4.353713116283091E-4-t33*t58*t75*2.846736678889049E-4+t33*t59*t74*8.255702643738717E-4-t57*t58*t59*1.639979074277157E-2+t33*t74*t75*8.255702643738717E-4-t57*t58*t75*1.639979074277157E-2-t57*t59*t74*5.915341110698354E-4-t58*t59*t73*6.001206764507874E-4-t57*t74*t75*5.915341110698354E-4-t58*t73*t75*6.001206764507874E-4-t59*t73*t74*1.635922650635785E-2-t73*t74*t75*1.635922650635785E-2-t9*t10*t33*t34*3.156734948660587E-2-t9*t10*t34*t59*1.113095897646181E-3-t9*t10*t34*t75*1.113095897646181E-3+t10*t33*t34*t59*2.896401859571867E-2+t10*t33*t34*t75*2.896401859571867E-2-t9*t33*t58*t59*5.585324801683264E-2-t10*t34*t57*t59*5.576841714485674E-2-t9*t33*t58*t75*5.585324801683264E-2+t9*t33*t59*t74*5.568788485803049E-2-t10*t34*t57*t75*5.576841714485674E-2+t10*t34*t59*t73*5.577259406372051E-2+t9*t33*t74*t75*5.568788485803049E-2+t10*t34*t73*t75*5.577259406372051E-2+t9*t10*t33*t34*t59*3.156734948660587E-2+t9*t10*t33*t34*t75*3.156734948660587E-2;
    c[11] = t109+t112-t9*t58*8.817574215578565E-4+t9*t59*7.546580486902988E-4-t13*t60*1.81E-2+t9*t74*1.270993728675576E-4+t9*t75*7.546580486902988E-4-t13*t76*1.81E-2-t33*t58*2.99216464372371E-3+t33*t59*3.029602411709927E-3-t33*t74*3.743776798621798E-5+t33*t75*3.029602411709927E-3+t57*t58*8.240645190217156E-3+t57*t59*1.557524060708425E-3+t57*t74*9.801830749074416E-3+t58*t73*8.240645190217156E-3+t57*t75*1.557524060708425E-3+t59*t73*1.557524060708425E-3+t73*t74*9.801830749074416E-3+t73*t75*1.557524060708425E-3-t9*t10*t34*2.959208971395655E-3+t10*t33*t34*1.009577116650697E-3-t10*t11*t57*6.574122182670539E-4-t10*t11*t73*6.574122182670539E-4+t10*t34*t57*2.297574320158763E-5-t11*t34*t57*1.549973690114125E-2+t10*t34*t73*2.297574320158763E-5-t11*t34*t73*1.549973690114125E-2+t9*t58*t59*1.270993728675576E-4+t9*t58*t75*1.270993728675576E-4-t9*t59*t74*8.817574215578565E-4-t33*t58*t59*3.743776798621798E-5-t35*t57*t58*2.875818595898751E-1-t9*t74*t75*8.817574215578565E-4-t33*t58*t75*3.743776798621798E-5-t33*t59*t74*2.99216464372371E-3-t35*t57*t74*2.875818595898751E-1-t35*t58*t73*2.875818595898751E-1+t57*t58*t59*9.801830749074416E-3-t33*t74*t75*2.99216464372371E-3-t35*t73*t74*2.875818595898751E-1+t57*t58*t75*9.801830749074416E-3+t57*t59*t74*8.240645190217156E-3+t58*t59*t73*9.801830749074416E-3+t57*t74*t75*8.240645190217156E-3+t58*t73*t75*9.801830749074416E-3+t59*t73*t74*8.240645190217156E-3+t73*t74*t75*8.240645190217156E-3+t9*t10*t34*t59*2.959208971395655E-3+t9*t10*t34*t75*2.959208971395655E-3-t10*t33*t34*t59*1.009577116650697E-3-t10*t33*t34*t75*1.009577116650697E-3-t10*t34*t57*t59*2.297574320158763E-5-t10*t34*t57*t75*2.297574320158763E-5-t10*t34*t59*t73*2.297574320158763E-5-t10*t34*t73*t75*2.297574320158763E-5;
    c[12] = t16*-3.370300000000001E-2+t40*2.47718E-1-t62*1.98970623391096E-4-t63*2.003180563211872E-5+t78*2.190024290232146E-4-t79*2.003180563211872E-5+t14*t15*6.725214316796237E-2+t16*t17*1.20005624647208E-1+t14*t38*1.586414893789752E-3-t15*t38*4.954563135453205E-1-t16*t41*6.370345907739961E-5-t17*t40*6.370345907739961E-5-t16*t65*3.370300000000001E-2-t40*t41*1.20005624647208E-1-t16*t81*3.3703E-2-t39*t62*4.365860704565494E-4+t40*t65*2.47718E-1-t39*t78*4.365860704565494E-4+t40*t81*2.47718E-1+t62*t63*2.190024290232146E-4-t64*t65*6.0E-2+t62*t79*2.190024290232146E-4-t63*t78*1.98970623391096E-4-t64*t81*6.0E-2-t65*t80*6.0E-2-t78*t79*1.98970623391096E-4-t80*t81*6.0E-2-t14*t38*t63*1.586414893789752E-3-t14*t38*t79*1.586414893789752E-3-6.0E-2;
    c[13] = t16*-2.47718E-1-t40*3.370300000000001E-2+t62*2.695145178952449E-5-t63*1.636252192216505E-3+t78*1.60930074042698E-3-t79*1.636252192216505E-3+t14*t15*4.941905177865888E-1+t16*t17*6.370345907739961E-5-t14*t38*4.170075425381788E-4+t15*t38*6.711175107535902E-2+t16*t41*1.20005624647208E-1+t17*t40*1.20005624647208E-1-t16*t65*2.47718E-1-t40*t41*6.370345907739961E-5-t16*t81*2.47718E-1-t39*t62*3.566153386244494E-2-t40*t65*3.370300000000001E-2-t39*t78*3.566153386244494E-2-t40*t81*3.3703E-2+t62*t63*1.60930074042698E-3+t62*t79*1.60930074042698E-3+t63*t78*2.695145178952449E-5+t78*t79*2.695145178952449E-5+t14*t38*t63*4.170075425381788E-4+t14*t38*t79*4.170075425381788E-4;
    c[14] = t62*2.299824701031743E-2+t63*1.170180756433687E-4+t78*2.288473491403919E-2+t79*1.170180756433687E-4-t14*t15*3.539606431708408E-2+t14*t38*2.84294570084937E-5-t15*t38*4.365115769377904E-3-t39*t62*4.987264424463382E-1-t39*t78*4.987264424463382E-1+t62*t63*2.288473491403919E-2-t64*t65*1.051397136175053E-2+t62*t79*2.288473491403919E-2+t63*t78*2.299824701031743E-2-t64*t81*1.051397136175053E-2-t65*t80*1.051397136175053E-2+t78*t79*2.299824701031743E-2-t80*t81*1.051397136175053E-2-t14*t38*t63*2.84294570084937E-5-t14*t38*t79*2.84294570084937E-5;
    c[15] = t4*-3.370300000000001E-2-t28*2.47718E-1-t50*1.98970623391096E-4-t51*2.003180563211872E-5+t66*2.190024290232146E-4-t67*2.003180563211872E-5+t2*t3*6.725214316796237E-2+t4*t5*1.20005624647208E-1-t2*t26*1.586414893789752E-3+t3*t26*4.954563135453205E-1+t4*t29*6.370345907739961E-5+t5*t28*6.370345907739961E-5-t4*t53*3.370300000000001E-2-t28*t29*1.20005624647208E-1-t4*t69*3.3703E-2-t27*t50*4.365860704565494E-4-t28*t53*2.47718E-1-t27*t66*4.365860704565494E-4-t28*t69*2.47718E-1+t50*t51*2.190024290232146E-4-t52*t53*6.0E-2+t50*t67*2.190024290232146E-4-t51*t66*1.98970623391096E-4-t52*t69*6.0E-2-t53*t68*6.0E-2-t66*t67*1.98970623391096E-4-t68*t69*6.0E-2+t2*t26*t51*1.586414893789752E-3+t2*t26*t67*1.586414893789752E-3-6.0E-2;
    c[16] = t4*2.47718E-1-t28*3.370300000000001E-2-t50*2.695145178952449E-5+t51*1.636252192216505E-3-t66*1.60930074042698E-3+t67*1.636252192216505E-3-t2*t3*4.941905177865888E-1-t4*t5*6.370345907739961E-5-t2*t26*4.170075425381788E-4+t3*t26*6.711175107535902E-2+t4*t29*1.20005624647208E-1+t5*t28*1.20005624647208E-1+t4*t53*2.47718E-1+t28*t29*6.370345907739961E-5+t4*t69*2.47718E-1+t27*t50*3.566153386244494E-2-t28*t53*3.370300000000001E-2+t27*t66*3.566153386244494E-2-t28*t69*3.3703E-2-t50*t51*1.60930074042698E-3-t50*t67*1.60930074042698E-3-t51*t66*2.695145178952449E-5-t66*t67*2.695145178952449E-5+t2*t26*t51*4.170075425381788E-4+t2*t26*t67*4.170075425381788E-4;
    c[17] = t50*2.299824701031743E-2+t51*1.170180756433687E-4+t66*2.288473491403919E-2+t67*1.170180756433687E-4-t2*t3*3.539606431708408E-2-t2*t26*2.84294570084937E-5+t3*t26*4.365115769377904E-3-t27*t50*4.987264424463382E-1-t27*t66*4.987264424463382E-1+t50*t51*2.288473491403919E-2-t52*t53*1.051397136175053E-2+t50*t67*2.288473491403919E-2+t51*t66*2.299824701031743E-2-t52*t69*1.051397136175053E-2-t53*t68*1.051397136175053E-2+t66*t67*2.299824701031743E-2-t68*t69*1.051397136175053E-2+t2*t26*t51*2.84294570084937E-5+t2*t26*t67*2.84294570084937E-5;

    fkPtr_->compute(0, 
                    contact_joint_id, 
                    q, 
                    nullptr, 
                    &stance_foot_endT, 
                    0);

    c.segment(18, 3) = fkPtr_->getTranslation() - stance_foot_T_des.p;
    c.segment(21, 3) = fkPtr_->getRPY() - stance_foot_T_des.getRPY();
}

void DigitDynamicsConstraints::get_J(const VecX& q) {
    assert(J.rows() == NUM_DEPENDENT_JOINTS);
    assert(J.cols() == modelPtr_->nv);

    double t2 = cos(q(9));
    double t3 = cos(q(10));
    double t4 = cos(q(11));
    double t5 = cos(q(12));
    double t6 = cos(q(13));
    double t7 = cos(q(14));
    double t8 = cos(q(15));
    double t9 = cos(q(16));
    double t10 = cos(q(17));
    double t11 = cos(q(18));
    double t12 = cos(q(19));
    double t13 = cos(q(20));
    double t14 = cos(q(24));
    double t15 = cos(q(25));
    double t16 = cos(q(26));
    double t17 = cos(q(27));
    double t18 = cos(q(28));
    double t19 = cos(q(29));
    double t20 = cos(q(30));
    double t21 = cos(q(31));
    double t22 = cos(q(32));
    double t23 = cos(q(33));
    double t24 = cos(q(34));
    double t25 = cos(q(35));
    double t26 = sin(q(9));
    double t27 = sin(q(10));
    double t28 = sin(q(11));
    double t29 = sin(q(12));
    double t30 = sin(q(13));
    double t31 = sin(q(14));
    double t32 = sin(q(15));
    double t33 = sin(q(16));
    double t34 = sin(q(17));
    double t35 = sin(q(18));
    double t36 = sin(q(19));
    double t37 = sin(q(20));
    double t38 = sin(q(24));
    double t39 = sin(q(25));
    double t40 = sin(q(26));
    double t41 = sin(q(27));
    double t42 = sin(q(28));
    double t43 = sin(q(29));
    double t44 = sin(q(30));
    double t45 = sin(q(31));
    double t46 = sin(q(32));
    double t47 = sin(q(33));
    double t48 = sin(q(34));
    double t49 = sin(q(35));
    double t50 = t2*t2;
    double t51 = t3*t3;
    double t52 = t5*t5;
    double t53 = t6*t6;
    double t54 = t7*t7;
    double t55 = t8*t8;
    double t56 = t9*t9;
    double t57 = t10*t10;
    double t58 = t11*t11;
    double t59 = t12*t12;
    double t60 = t13*t13;
    double t61 = t14*t14;
    double t62 = t15*t15;
    double t63 = t17*t17;
    double t64 = t26*t26;
    double t65 = t27*t27;
    double t66 = t29*t29;
    double t67 = t30*t30;
    double t68 = t31*t31;
    double t69 = t32*t32;
    double t70 = t33*t33;
    double t71 = t34*t34;
    double t72 = t35*t35;
    double t73 = t36*t36;
    double t74 = t37*t37;
    double t75 = t38*t38;
    double t76 = t39*t39;
    double t77 = t41*t41;
    double t78 = t24*1.696216709330505E-2;
    double t79 = t25*-9.551E-3;
    double t80 = t25*9.551E-3;
    double t81 = t48*-1.696216709330505E-2;
    double t82 = t48*1.696216709330505E-2;
    double t83 = t24*5.143951577823025E-2;
    double t84 = t48*5.143951577823025E-2;
    double t89 = t4*t5*6.370345907739961E-5;
    double t90 = t16*t17*-6.370345907739961E-5;
    double t91 = t16*t17*6.370345907739961E-5;
    double t92 = t4*t29*6.370345907739961E-5;
    double t93 = t5*t28*6.370345907739961E-5;
    double t94 = t16*t41*-6.370345907739961E-5;
    double t95 = t16*t41*6.370345907739961E-5;
    double t96 = t17*t40*-6.370345907739961E-5;
    double t97 = t17*t40*6.370345907739961E-5;
    double t98 = t28*t29*-6.370345907739961E-5;
    double t99 = t28*t29*6.370345907739961E-5;
    double t100 = t40*t41*6.370345907739961E-5;
    double t101 = t12*t36*3.517668160995551E-1;
    double t102 = t4*t5*1.20005624647208E-1;
    double t103 = t16*t17*1.20005624647208E-1;
    double t104 = t12*t36*-2.105329121329532E-1;
    double t105 = t12*t36*2.105329121329532E-1;
    double t106 = t4*t29*-1.20005624647208E-1;
    double t107 = t4*t29*1.20005624647208E-1;
    double t108 = t5*t28*-1.20005624647208E-1;
    double t109 = t5*t28*1.20005624647208E-1;
    double t110 = t16*t41*-1.20005624647208E-1;
    double t111 = t16*t41*1.20005624647208E-1;
    double t112 = t17*t40*-1.20005624647208E-1;
    double t113 = t17*t40*1.20005624647208E-1;
    double t114 = t28*t29*-1.20005624647208E-1;
    double t115 = t28*t29*1.20005624647208E-1;
    double t116 = t40*t41*-1.20005624647208E-1;
    double t117 = t40*t41*1.20005624647208E-1;
    double t85 = t73*1.758834080497775E-1;
    double t86 = t59*-1.052664560664766E-1;
    double t87 = t59*1.052664560664766E-1;
    double t88 = t73*1.052664560664766E-1;
    double t118 = t59*-1.758834080497775E-1;
    double t119 = t59*1.758834080497775E-1;
    double t120 = t12*t60*1.696216709330505E-2;
    double t121 = t13*t59*9.551E-3;
    double t122 = t12*t74*1.696216709330505E-2;
    double t123 = t36*t60*-1.696216709330505E-2;
    double t124 = t36*t60*1.696216709330505E-2;
    double t125 = t13*t73*9.551E-3;
    double t126 = t36*t74*-1.696216709330505E-2;
    double t127 = t36*t74*1.696216709330505E-2;
    double t128 = t12*t60*-5.143951577823025E-2;
    double t129 = t12*t60*5.143951577823025E-2;
    double t130 = t12*t74*-5.143951577823025E-2;
    double t131 = t12*t74*5.143951577823025E-2;
    double t132 = t36*t60*-5.143951577823025E-2;
    double t133 = t36*t60*5.143951577823025E-2;
    double t134 = t36*t74*-5.143951577823025E-2;
    double t135 = t36*t74*5.143951577823025E-2;
    double t138 = t60*t73*-1.758834080497775E-1;
    double t141 = t73*t74*-1.758834080497775E-1;
    double t144 = t60*t73*-1.052664560664766E-1;
    double t146 = t73*t74*-1.052664560664766E-1;
    double t148 = t12*t36*t60*-3.517668160995551E-1;
    double t149 = t60*t101;
    double t150 = t12*t36*t74*-3.517668160995551E-1;
    double t151 = t74*t101;
    double t152 = t60*t105;
    double t153 = t74*t105;
    double t136 = t60*t119;
    double t137 = t74*t119;
    double t139 = t60*t85;
    double t140 = t60*t87;
    double t142 = t74*t85;
    double t143 = t74*t87;
    double t145 = t60*t88;
    double t147 = t74*t88;

    // assume J(NUM_DEPENDENT_JOINTS, modelPtr_->nv) has been allocated outside
    J.setZero();
    J(0, 28) = t18*5.699060997402858E-2-t42*1.034589188110661E-3-t18*t44*1.26575449899492E-2+t42*t44*2.246221860400801E-3+t18*t19*t20*3.020283789547073E-3-t18*t20*t43*3.397508858570615E-1-t19*t20*t42*3.399783924207052E-1-t20*t42*t43*3.105990081579729E-3;
    J(0, 29) = t18*t19*t20*3.105990081579729E-3-t18*t20*t43*3.399783924207052E-1-t19*t20*t42*3.397508858570615E-1-t20*t42*t43*3.020283789547073E-3;
    J(0, 30) = t18*t20*-2.246221860400801E-3-t20*t42*1.26575449899492E-2-t18*t19*t44*3.399783924207052E-1-t18*t43*t44*3.105990081579729E-3-t19*t42*t44*3.020283789547073E-3+t42*t43*t44*3.397508858570615E-1;
    J(0, 34) = t81+t83+t24*t25*2.991020934719675E-3+t24*t49*5.605619802270151E-3+t25*t48*9.070578524442013E-3+t48*t49*1.69996184260823E-2;
    J(0, 35) = t24*t25*-1.69996184260823E-2+t24*t49*9.070578524442013E-3+t25*t48*5.605619802270151E-3-t48*t49*2.991020934719675E-3;
    J(1, 28) = t18*-1.034589188110661E-3-t42*5.699060997402858E-2+t18*t44*2.246221860400801E-3+t42*t44*1.26575449899492E-2-t18*t19*t20*3.399783924207052E-1-t18*t20*t43*3.105990081579729E-3-t19*t20*t42*3.020283789547073E-3+t20*t42*t43*3.397508858570615E-1;
    J(1, 29) = t18*t19*t20*-3.397508858570615E-1-t18*t20*t43*3.020283789547073E-3-t19*t20*t42*3.105990081579729E-3+t20*t42*t43*3.399783924207052E-1;
    J(1, 30) = t18*t20*-1.26575449899492E-2+t20*t42*2.246221860400801E-3-t18*t19*t44*3.020283789547073E-3+t18*t43*t44*3.397508858570615E-1+t19*t42*t44*3.399783924207052E-1+t42*t43*t44*3.105990081579729E-3;
    J(1, 34) = t78+t84-t24*t25*9.070578524442013E-3-t24*t49*1.69996184260823E-2+t25*t48*2.991020934719675E-3+t48*t49*5.605619802270151E-3;
    J(1, 35) = t24*t25*-5.605619802270151E-3+t24*t49*2.991020934719675E-3-t25*t48*1.69996184260823E-2+t48*t49*9.070578524442013E-3;
    J(2, 29) = t19*t20*-1.263678697118524E-2-t20*t43*2.360206106172353E-3;
    J(2, 30) = t20*3.39756885202024E-1-t19*t44*2.360206106172353E-3+t43*t44*1.263678697118524E-2;
    J(2, 35) = t49*-1.79E-2+t79;
    J(3, 31) = t21*-5.699060997402856E-2-t45*1.034589188110661E-3-t21*t47*1.549180159908108E-2+t45*t47*8.234792310892687E-4-t21*t22*t23*3.104080344633556E-3-t21*t23*t46*2.87566285868891E-1-t22*t23*t45*2.879825211612492E-1+t23*t45*t46*3.064210757541298E-3;
    J(3, 32) = t21*t22*t23*-3.064210757541298E-3-t21*t23*t46*2.879825211612492E-1-t22*t23*t45*2.87566285868891E-1+t23*t45*t46*3.104080344633556E-3;
    J(3, 33) = t21*t23*-8.234792310892687E-4-t23*t45*1.549180159908108E-2-t21*t22*t47*2.879825211612492E-1+t21*t46*t47*3.064210757541298E-3+t22*t45*t47*3.104080344633556E-3+t45*t46*t47*2.87566285868891E-1;
    J(3, 34) = t81+t83+t24*t25*2.991020934719677E-3-t24*t49*5.668252425759201E-3+t25*t48*9.070578524442012E-3-t48*t49*1.718955829676478E-2;
    J(3, 35) = t24*t25*1.718955829676478E-2+t24*t49*9.070578524442012E-3-t25*t48*5.668252425759201E-3-t48*t49*2.991020934719677E-3;
    J(4, 31) = t21*1.034589188110661E-3-t45*5.699060997402856E-2-t21*t47*8.234792310892687E-4-t45*t47*1.549180159908108E-2+t21*t22*t23*2.879825211612492E-1-t21*t23*t46*3.064210757541298E-3-t22*t23*t45*3.104080344633556E-3-t23*t45*t46*2.87566285868891E-1;
    J(4, 32) = t21*t22*t23*2.87566285868891E-1-t21*t23*t46*3.104080344633556E-3-t22*t23*t45*3.064210757541298E-3-t23*t45*t46*2.879825211612492E-1;
    J(4, 33) = t21*t23*1.549180159908108E-2-t23*t45*8.234792310892687E-4-t21*t22*t47*3.104080344633556E-3-t21*t46*t47*2.87566285868891E-1-t22*t45*t47*2.879825211612492E-1+t45*t46*t47*3.064210757541298E-3;
    J(4, 34) = t78+t84-t24*t25*9.070578524442012E-3+t24*t49*1.718955829676478E-2+t25*t48*2.991020934719677E-3-t48*t49*5.668252425759201E-3;
    J(4, 35) = t24*t25*5.668252425759201E-3+t24*t49*2.991020934719677E-3+t25*t48*1.718955829676478E-2+t48*t49*9.070578524442012E-3;
    J(5, 32) = t22*t23*1.549973690114125E-2+t23*t46*6.574122182670539E-4;
    J(5, 33) = t23*-2.875818595898751E-1+t22*t47*6.574122182670539E-4-t46*t47*1.549973690114125E-2;
    J(5, 35) = t49*1.81E-2+t79;
    J(6, 13) = t6*t54*-4.104331657223103E-4-t6*t55*2.809111102928818E-2-t6*t68*2.848906577901809E-2-t6*t69*2.809111102928818E-2-t30*t54*7.74045457557333E-4-t30*t55*4.455658896021977E-4+t30*t68*1.8502215904887E-4-t30*t69*4.455658896021977E-4+t53*t54*1.643509328293678E-2+t53*t55*3.732977901401361E-5-t53*t68*1.64724230619508E-2-t54*t67*1.643509328293678E-2+t53*t69*3.732977901401361E-5-t55*t67*3.732977901401361E-5+t67*t68*1.64724230619508E-2-t67*t69*3.732977901401361E-5-t6*t7*t8*3.020283789547073E-3-t6*t7*t31*9.515128902424641E-4-t6*t8*t31*3.397508858570615E-1-t7*t8*t30*3.399783924207052E-1-t7*t30*t31*2.809726196867663E-2+t8*t30*t31*3.105990081579729E-3-t6*t30*t54*5.952939355245877E-2-t6*t30*t55*6.247020354403988E-5-t7*t31*t53*5.956062305521773E-2+t6*t32*t54*1.26575449899492E-2+t6*t30*t68*5.959186375600283E-2-t6*t30*t69*6.247020354403988E-5+t7*t31*t67*5.956062305521773E-2+t6*t32*t68*1.26575449899492E-2-t6*t54*t55*2.848906577901809E-2+t30*t32*t54*2.246221860400801E-3-t6*t54*t69*2.848906577901809E-2-t6*t55*t68*4.104331657223103E-4+t30*t32*t68*2.246221860400801E-3+t30*t54*t55*1.8502215904887E-4-t6*t68*t69*4.104331657223103E-4+t30*t54*t69*1.8502215904887E-4-t30*t55*t68*7.74045457557333E-4-t53*t54*t55*1.64724230619508E-2-t30*t68*t69*7.74045457557333E-4-t53*t54*t69*1.64724230619508E-2+t53*t55*t68*1.643509328293678E-2+t54*t55*t67*1.64724230619508E-2+t53*t68*t69*1.643509328293678E-2+t54*t67*t69*1.64724230619508E-2-t55*t67*t68*1.643509328293678E-2-t67*t68*t69*1.643509328293678E-2-t6*t7*t30*t31*6.581499972336942E-2+t6*t7*t31*t55*9.515128902424641E-4+t6*t7*t31*t69*9.515128902424641E-4+t7*t30*t31*t55*2.809726196867663E-2+t7*t30*t31*t69*2.809726196867663E-2+t6*t30*t54*t55*5.959186375600283E-2+t7*t31*t53*t55*5.956062305521773E-2+t6*t30*t54*t69*5.959186375600283E-2-t6*t30*t55*t68*5.952939355245877E-2+t7*t31*t53*t69*5.956062305521773E-2-t7*t31*t55*t67*5.956062305521773E-2-t6*t30*t68*t69*5.952939355245877E-2-t7*t31*t67*t69*5.956062305521773E-2+t6*t7*t30*t31*t55*6.581499972336942E-2+t6*t7*t30*t31*t69*6.581499972336942E-2;
    J(6, 14) = t54*-1.522231734027378E-5+t68*1.522231734027378E-5-t7*t31*7.865874492559456E-5+t6*t54*2.809726196867663E-2-t6*t68*2.809726196867663E-2-t30*t54*9.515128902424641E-4+t30*t68*9.515128902424641E-4+t53*t54*1.646136108951249E-2+t54*t55*1.522231734027378E-5-t53*t68*1.646136108951249E-2-t54*t67*1.644613877217222E-2+t54*t69*1.522231734027378E-5-t55*t68*1.522231734027378E-5+t67*t68*1.644613877217222E-2-t68*t69*1.522231734027378E-5-t6*t7*t8*3.105990081579729E-3-t6*t7*t31*1.918135233212406E-3-t6*t8*t31*3.399783924207052E-1-t7*t8*t30*3.397508858570615E-1-t7*t30*t31*5.615726522659155E-2+t8*t30*t31*3.020283789547073E-3-t6*t30*t54*5.956062305521773E-2-t7*t31*t53*5.952129928176799E-2+t7*t31*t55*7.865874492559456E-5+t6*t30*t68*5.956062305521773E-2+t7*t31*t67*5.959995802669361E-2+t7*t31*t69*7.865874492559456E-5-t6*t54*t55*2.809726196867663E-2-t6*t54*t69*2.809726196867663E-2+t6*t55*t68*2.809726196867663E-2+t30*t54*t55*9.515128902424641E-4+t6*t68*t69*2.809726196867663E-2+t30*t54*t69*9.515128902424641E-4-t30*t55*t68*9.515128902424641E-4-t53*t54*t55*1.646136108951249E-2-t30*t68*t69*9.515128902424641E-4-t53*t54*t69*1.646136108951249E-2+t53*t55*t68*1.646136108951249E-2+t54*t55*t67*1.644613877217222E-2+t53*t68*t69*1.646136108951249E-2+t54*t67*t69*1.644613877217222E-2-t55*t67*t68*1.644613877217222E-2-t67*t68*t69*1.644613877217222E-2-t6*t7*t30*t31*6.581503268977515E-2+t6*t7*t31*t55*1.918135233212406E-3+t6*t7*t31*t69*1.918135233212406E-3+t7*t30*t31*t55*5.615726522659155E-2+t7*t30*t31*t69*5.615726522659155E-2+t6*t30*t54*t55*5.956062305521773E-2+t7*t31*t53*t55*5.952129928176799E-2+t6*t30*t54*t69*5.956062305521773E-2-t6*t30*t55*t68*5.956062305521773E-2+t7*t31*t53*t69*5.952129928176799E-2-t7*t31*t55*t67*5.959995802669361E-2-t6*t30*t68*t69*5.956062305521773E-2-t7*t31*t67*t69*5.959995802669361E-2+t6*t7*t30*t31*t55*6.581503268977515E-2+t6*t7*t30*t31*t69*6.581503268977515E-2;
    J(6, 15) = t6*t7*t32*-3.399783924207052E-1-t6*t8*t54*2.246221860400801E-3+t6*t31*t32*3.105990081579729E-3+t7*t30*t32*3.020283789547073E-3-t6*t8*t68*2.246221860400801E-3+t8*t30*t54*1.26575449899492E-2+t30*t31*t32*3.397508858570615E-1+t8*t30*t68*1.26575449899492E-2;
    J(6, 19) = t86+t88+t101+t123+t126+t128+t130+t140+t143+t144+t146+t148+t150-t12*t13*2.991020934719675E-3+t12*t37*5.605619802270151E-3+t13*t36*9.070578524442013E-3-t36*t37*1.69996184260823E-2;
    J(6, 20) = t12*t13*1.69996184260823E-2+t12*t37*9.070578524442013E-3+t13*t36*5.605619802270151E-3+t36*t37*2.991020934719675E-3;
    J(7, 13) = t6*t54*-7.74045457557333E-4-t6*t55*4.455658896021977E-4+t6*t68*1.8502215904887E-4-t6*t69*4.455658896021977E-4+t30*t54*4.104331657223105E-4+t30*t55*2.809111102928818E-2+t30*t68*2.848906577901809E-2+t30*t69*2.809111102928818E-2-t53*t54*2.976469677622939E-2-t53*t55*3.123510177201994E-5+t53*t68*2.979593187800142E-2+t54*t67*2.976469677622939E-2-t53*t69*3.123510177201994E-5+t55*t67*3.123510177201994E-5-t67*t68*2.979593187800142E-2+t67*t69*3.123510177201994E-5-t6*t7*t8*3.399783924207052E-1-t6*t7*t31*2.809726196867663E-2+t6*t8*t31*3.105990081579729E-3+t7*t8*t30*3.020283789547073E-3+t7*t30*t31*9.51512890242464E-4+t8*t30*t31*3.397508858570615E-1-t6*t30*t54*3.287018656587355E-2-t6*t30*t55*7.465955802802721E-5-t7*t31*t53*3.290749986168471E-2+t6*t32*t54*2.246221860400801E-3+t6*t30*t68*3.294484612390159E-2-t6*t30*t69*7.465955802802721E-5+t7*t31*t67*3.290749986168471E-2+t6*t32*t68*2.246221860400801E-3+t6*t54*t55*1.8502215904887E-4-t30*t32*t54*1.26575449899492E-2+t6*t54*t69*1.8502215904887E-4-t6*t55*t68*7.74045457557333E-4-t30*t32*t68*1.26575449899492E-2+t30*t54*t55*2.848906577901809E-2-t6*t68*t69*7.74045457557333E-4+t30*t54*t69*2.848906577901809E-2+t30*t55*t68*4.104331657223105E-4+t53*t54*t55*2.979593187800142E-2+t30*t68*t69*4.104331657223105E-4+t53*t54*t69*2.979593187800142E-2-t53*t55*t68*2.976469677622939E-2-t54*t55*t67*2.979593187800142E-2-t53*t68*t69*2.976469677622939E-2-t54*t67*t69*2.979593187800142E-2+t55*t67*t68*2.976469677622939E-2+t67*t68*t69*2.976469677622939E-2+t6*t7*t30*t31*1.191212461104355E-1+t6*t7*t31*t55*2.809726196867663E-2+t6*t7*t31*t69*2.809726196867663E-2-t7*t30*t31*t55*9.51512890242464E-4-t7*t30*t31*t69*9.51512890242464E-4+t6*t30*t54*t55*3.294484612390159E-2+t7*t31*t53*t55*3.290749986168471E-2+t6*t30*t54*t69*3.294484612390159E-2-t6*t30*t55*t68*3.287018656587355E-2+t7*t31*t53*t69*3.290749986168471E-2-t7*t31*t55*t67*3.290749986168471E-2-t6*t30*t68*t69*3.287018656587355E-2-t7*t31*t67*t69*3.290749986168471E-2-t6*t7*t30*t31*t55*1.191212461104355E-1-t6*t7*t30*t31*t69*1.191212461104355E-1;
    J(7, 14) = t54*8.772182874056074E-6-t68*8.772182874056074E-6+t7*t31*4.532876826220704E-5-t6*t54*9.51512890242464E-4+t6*t68*9.51512890242464E-4-t30*t54*2.809726196867663E-2+t30*t68*2.809726196867663E-2-t53*t54*2.978469761904589E-2-t54*t55*8.772182874056074E-6+t53*t68*2.978469761904589E-2+t54*t67*2.977592543617184E-2-t54*t69*8.772182874056074E-6+t55*t68*8.772182874056074E-6-t67*t68*2.977592543617184E-2+t68*t69*8.772182874056074E-6-t6*t7*t8*3.397508858570615E-1-t6*t7*t31*5.615726522659155E-2+t6*t8*t31*3.020283789547073E-3+t7*t8*t30*3.105990081579729E-3+t7*t30*t31*1.918135233212406E-3+t8*t30*t31*3.399783924207052E-1-t6*t30*t54*3.290749986168471E-2-t7*t31*t53*3.293018072901868E-2-t7*t31*t55*4.532876826220704E-5+t6*t30*t68*3.290749986168471E-2+t7*t31*t67*3.288485196075646E-2-t7*t31*t69*4.532876826220704E-5+t6*t54*t55*9.51512890242464E-4+t6*t54*t69*9.51512890242464E-4-t6*t55*t68*9.51512890242464E-4+t30*t54*t55*2.809726196867663E-2-t6*t68*t69*9.51512890242464E-4+t30*t54*t69*2.809726196867663E-2-t30*t55*t68*2.809726196867663E-2+t53*t54*t55*2.978469761904589E-2-t30*t68*t69*2.809726196867663E-2+t53*t54*t69*2.978469761904589E-2-t53*t55*t68*2.978469761904589E-2-t54*t55*t67*2.977592543617184E-2-t53*t68*t69*2.978469761904589E-2-t54*t67*t69*2.977592543617184E-2+t55*t67*t68*2.977592543617184E-2+t67*t68*t69*2.977592543617184E-2+t6*t7*t30*t31*1.191212573084616E-1+t6*t7*t31*t55*5.615726522659155E-2+t6*t7*t31*t69*5.615726522659155E-2-t7*t30*t31*t55*1.918135233212406E-3-t7*t30*t31*t69*1.918135233212406E-3+t6*t30*t54*t55*3.290749986168471E-2+t7*t31*t53*t55*3.293018072901868E-2+t6*t30*t54*t69*3.290749986168471E-2-t6*t30*t55*t68*3.290749986168471E-2+t7*t31*t53*t69*3.293018072901868E-2-t7*t31*t55*t67*3.288485196075646E-2-t6*t30*t68*t69*3.290749986168471E-2-t7*t31*t67*t69*3.288485196075646E-2-t6*t7*t30*t31*t55*1.191212573084616E-1-t6*t7*t30*t31*t69*1.191212573084616E-1;
    J(7, 15) = t6*t7*t32*3.020283789547073E-3+t6*t8*t54*1.26575449899492E-2+t6*t31*t32*3.397508858570615E-1+t7*t30*t32*3.399783924207052E-1+t6*t8*t68*1.26575449899492E-2+t8*t30*t54*2.246221860400801E-3-t30*t31*t32*3.105990081579729E-3+t8*t30*t68*2.246221860400801E-3;
    J(7, 19) = t85+t104+t118+t120+t122+t132+t134+t136+t137+t138+t141+t152+t153-t12*t13*9.070578524442013E-3+t12*t37*1.69996184260823E-2-t13*t36*2.991020934719675E-3+t36*t37*5.605619802270151E-3;
    J(7, 20) = t12*t13*-5.605619802270151E-3-t12*t37*2.991020934719675E-3+t13*t36*1.69996184260823E-2+t36*t37*9.070578524442013E-3;
    J(8, 13) = t6*t54*1.101395785003593E-3-t6*t55*9.852121016062406E-4-t6*t68*1.161836833973452E-4-t6*t69*9.852121016062406E-4+t30*t54*6.213602549397351E-4-t30*t55*8.271781369638137E-4+t30*t68*2.058178820240719E-4-t30*t69*8.271781369638137E-4+t6*t7*t31*4.163488242625716E-4-t7*t30*t31*1.218023269812134E-3+t6*t30*t54*1.040834085586084E-17-t6*t30*t68*1.040834085586084E-17-t6*t54*t55*1.161836833973452E-4-t6*t54*t69*1.161836833973452E-4+t6*t55*t68*1.101395785003593E-3+t30*t54*t55*2.058178820240719E-4+t6*t68*t69*1.101395785003593E-3+t30*t54*t69*2.058178820240719E-4+t30*t55*t68*6.213602549397351E-4+t30*t68*t69*6.213602549397351E-4-t6*t7*t30*t31*6.396792817664476E-18-t6*t7*t31*t55*4.163488242625716E-4-t6*t7*t31*t69*4.163488242625716E-4+t7*t30*t31*t55*1.218023269812134E-3+t7*t30*t31*t69*1.218023269812134E-3-t6*t30*t54*t55*1.040834085586084E-17-t6*t30*t54*t69*1.040834085586084E-17+t6*t30*t55*t68*1.040834085586084E-17+t6*t30*t68*t69*1.040834085586084E-17+t6*t7*t30*t31*t55*6.396792817664476E-18+t6*t7*t30*t31*t69*6.396792817664476E-18;
    J(8, 14) = t6*t54*1.218023269812134E-3-t6*t68*1.218023269812134E-3+t30*t54*4.163488242625716E-4-t30*t68*4.163488242625716E-4+t53*t54*2.212067055639549E-4-t53*t68*2.212067055639549E-4+t54*t67*2.212067055639517E-4-t67*t68*2.212067055639517E-4+t6*t7*t31*8.310847458313263E-4+t7*t8*t53*1.263678697118524E-2-t7*t30*t31*2.435158936801877E-3+t7*t8*t67*1.263678697118524E-2+t7*t31*t53*2.079441461225871E-3-t8*t31*t53*2.360206106172353E-3+t7*t31*t67*2.07944146122585E-3-t8*t31*t67*2.360206106172353E-3-t6*t54*t55*1.218023269812134E-3-t6*t54*t69*1.218023269812134E-3+t6*t55*t68*1.218023269812134E-3-t30*t54*t55*4.163488242625716E-4+t6*t68*t69*1.218023269812134E-3-t30*t54*t69*4.163488242625716E-4+t30*t55*t68*4.163488242625716E-4-t53*t54*t55*2.212067055639549E-4+t30*t68*t69*4.163488242625716E-4-t53*t54*t69*2.212067055639549E-4+t53*t55*t68*2.212067055639549E-4-t54*t55*t67*2.212067055639517E-4+t53*t68*t69*2.212067055639549E-4-t54*t67*t69*2.212067055639517E-4+t55*t67*t68*2.212067055639517E-4+t67*t68*t69*2.212067055639517E-4-t6*t7*t31*t55*8.310847458313263E-4-t6*t7*t31*t69*8.310847458313263E-4+t7*t30*t31*t55*2.435158936801877E-3+t7*t30*t31*t69*2.435158936801877E-3-t7*t31*t53*t55*2.079441461225871E-3-t7*t31*t53*t69*2.079441461225871E-3-t7*t31*t55*t67*2.07944146122585E-3-t7*t31*t67*t69*2.07944146122585E-3;
    J(8, 15) = t7*t32*t53*-2.360206106172353E-3-t7*t32*t67*2.360206106172353E-3+t8*t53*t54*3.39756885202024E-1-t31*t32*t53*1.263678697118524E-2+t8*t53*t68*3.39756885202024E-1+t8*t54*t67*3.39756885202024E-1-t31*t32*t67*1.263678697118524E-2+t8*t67*t68*3.39756885202024E-1;
    J(8, 20) = t121+t125-t37*t59*1.79E-2-t37*t73*1.79E-2;
    J(9, 16) = t9*t57*-4.353713116283091E-4+t9*t58*2.893932048270239E-2+t9*t71*2.848666080295449E-2+t9*t72*2.893932048270239E-2-t33*t57*8.255702643738717E-4-t33*t58*4.936925916256938E-4+t33*t71*2.846736678889049E-4-t33*t72*4.936925916256938E-4-t56*t57*1.576769239528801E-2-t56*t58*3.197767103277491E-5+t56*t71*1.579967006632079E-2+t57*t70*1.576769239528801E-2-t56*t72*3.197767103277491E-5+t58*t70*3.197767103277491E-5-t70*t71*1.579967006632079E-2+t70*t72*3.197767103277491E-5+t9*t10*t11*3.104080344633556E-3-t9*t10*t34*1.113095897646181E-3-t9*t11*t34*2.87566285868891E-1-t10*t11*t33*2.879825211612492E-1+t10*t33*t34*2.896401859571867E-2-t11*t33*t34*3.064210757541298E-3-t9*t33*t57*1.11375769716061E-1-t9*t33*t58*3.307263176043031E-4-t10*t34*t56*1.115410112085772E-1+t9*t35*t57*1.549180159908108E-2+t9*t33*t71*1.117064960336653E-1-t9*t33*t72*3.307263176043031E-4+t10*t34*t70*1.115410112085772E-1+t9*t35*t71*1.549180159908108E-2+t9*t57*t58*2.848666080295449E-2+t33*t35*t57*8.234792310892687E-4+t9*t57*t72*2.848666080295449E-2-t9*t58*t71*4.353713116283091E-4+t33*t35*t71*8.234792310892687E-4+t33*t57*t58*2.846736678889049E-4-t9*t71*t72*4.353713116283091E-4+t33*t57*t72*2.846736678889049E-4-t33*t58*t71*8.255702643738717E-4+t56*t57*t58*1.579967006632079E-2-t33*t71*t72*8.255702643738717E-4+t56*t57*t72*1.579967006632079E-2-t56*t58*t71*1.576769239528801E-2-t57*t58*t70*1.579967006632079E-2-t56*t71*t72*1.576769239528801E-2-t57*t70*t72*1.579967006632079E-2+t58*t70*t71*1.576769239528801E-2+t70*t71*t72*1.576769239528801E-2+t9*t10*t33*t34*6.313469897321174E-2+t9*t10*t34*t58*1.113095897646181E-3+t9*t10*t34*t72*1.113095897646181E-3-t10*t33*t34*t58*2.896401859571867E-2-t10*t33*t34*t72*2.896401859571867E-2+t9*t33*t57*t58*1.117064960336653E-1+t10*t34*t56*t58*1.115410112085772E-1+t9*t33*t57*t72*1.117064960336653E-1-t9*t33*t58*t71*1.11375769716061E-1+t10*t34*t56*t72*1.115410112085772E-1-t10*t34*t58*t70*1.115410112085772E-1-t9*t33*t71*t72*1.11375769716061E-1-t10*t34*t70*t72*1.115410112085772E-1-t9*t10*t33*t34*t58*6.313469897321174E-2-t9*t10*t33*t34*t72*6.313469897321174E-2;
    J(9, 17) = t57*-1.363641158467861E-5+t71*1.363641158467861E-5-t10*t34*3.209258234829029E-4-t9*t57*2.896401859571867E-2+t9*t71*2.896401859571867E-2-t33*t57*1.113095897646181E-3+t33*t71*1.113095897646181E-3-t56*t57*1.577685653751059E-2+t57*t58*1.363641158467861E-5+t56*t71*1.577685653751059E-2+t57*t70*1.579049294909527E-2+t57*t72*1.363641158467861E-5-t58*t71*1.363641158467861E-5-t70*t71*1.579049294909527E-2-t71*t72*1.363641158467861E-5+t9*t10*t11*3.064210757541298E-3-t9*t10*t34*2.220487864525553E-3-t9*t11*t34*2.879825211612492E-1-t10*t11*t33*2.87566285868891E-1+t10*t33*t34*5.784406422916559E-2-t11*t33*t34*3.104080344633556E-3-t9*t33*t57*1.115410112085772E-1-t10*t34*t56*1.113806699631217E-1+t10*t34*t58*3.209258234829029E-4+t9*t33*t71*1.115410112085772E-1+t10*t34*t70*1.117015957866046E-1+t10*t34*t72*3.209258234829029E-4+t9*t57*t58*2.896401859571867E-2+t9*t57*t72*2.896401859571867E-2-t9*t58*t71*2.896401859571867E-2+t33*t57*t58*1.113095897646181E-3-t9*t71*t72*2.896401859571867E-2+t33*t57*t72*1.113095897646181E-3-t33*t58*t71*1.113095897646181E-3+t56*t57*t58*1.577685653751059E-2-t33*t71*t72*1.113095897646181E-3+t56*t57*t72*1.577685653751059E-2-t56*t58*t71*1.577685653751059E-2-t57*t58*t70*1.579049294909527E-2-t56*t71*t72*1.577685653751059E-2-t57*t70*t72*1.579049294909527E-2+t58*t70*t71*1.579049294909527E-2+t70*t71*t72*1.579049294909527E-2+t9*t10*t33*t34*6.31347249232176E-2+t9*t10*t34*t58*2.220487864525553E-3+t9*t10*t34*t72*2.220487864525553E-3-t10*t33*t34*t58*5.784406422916559E-2-t10*t33*t34*t72*5.784406422916559E-2+t9*t33*t57*t58*1.115410112085772E-1+t10*t34*t56*t58*1.113806699631217E-1+t9*t33*t57*t72*1.115410112085772E-1-t9*t33*t58*t71*1.115410112085772E-1+t10*t34*t56*t72*1.113806699631217E-1-t10*t34*t58*t70*1.117015957866046E-1-t9*t33*t71*t72*1.115410112085772E-1-t10*t34*t70*t72*1.117015957866046E-1-t9*t10*t33*t34*t58*6.31347249232176E-2-t9*t10*t33*t34*t72*6.31347249232176E-2;
    J(9, 18) = t9*t10*t35*-2.879825211612492E-1-t9*t11*t57*8.234792310892687E-4-t9*t34*t35*3.064210757541298E-3-t10*t33*t35*3.104080344633556E-3-t9*t11*t71*8.234792310892687E-4+t11*t33*t57*1.549180159908108E-2+t33*t34*t35*2.87566285868891E-1+t11*t33*t71*1.549180159908108E-2;
    J(9, 19) = t86+t88+t101+t123+t126+t128+t130+t140+t143+t144+t146+t148+t150-t12*t13*2.991020934719677E-3-t12*t37*5.668252425759201E-3+t13*t36*9.070578524442012E-3+t36*t37*1.718955829676478E-2;
    J(9, 20) = t12*t13*-1.718955829676478E-2+t12*t37*9.070578524442012E-3-t13*t36*5.668252425759201E-3+t36*t37*2.991020934719677E-3;
    J(10, 16) = t9*t57*8.255702643738717E-4+t9*t58*4.936925916256938E-4-t9*t71*2.846736678889049E-4+t9*t72*4.936925916256938E-4-t33*t57*4.353713116283091E-4+t33*t58*2.893932048270239E-2+t33*t71*2.848666080295449E-2+t33*t72*2.893932048270239E-2+t56*t57*5.568788485803049E-2+t56*t58*1.653631588021515E-4-t56*t71*5.585324801683264E-2-t57*t70*5.568788485803049E-2+t56*t72*1.653631588021515E-4-t58*t70*1.653631588021515E-4+t70*t71*5.585324801683264E-2-t70*t72*1.653631588021515E-4+t9*t10*t11*2.879825211612492E-1-t9*t10*t34*2.896401859571867E-2+t9*t11*t34*3.064210757541298E-3+t10*t11*t33*3.104080344633556E-3-t10*t33*t34*1.113095897646181E-3-t11*t33*t34*2.87566285868891E-1-t9*t33*t57*3.153538479057603E-2-t9*t33*t58*6.395534206554983E-5-t10*t34*t56*3.156734948660587E-2-t9*t35*t57*8.234792310892687E-4+t9*t33*t71*3.159934013264157E-2-t9*t33*t72*6.395534206554983E-5+t10*t34*t70*3.156734948660587E-2-t9*t35*t71*8.234792310892687E-4-t9*t57*t58*2.846736678889049E-4+t33*t35*t57*1.549180159908108E-2-t9*t57*t72*2.846736678889049E-4+t9*t58*t71*8.255702643738717E-4+t33*t35*t71*1.549180159908108E-2+t33*t57*t58*2.848666080295449E-2+t9*t71*t72*8.255702643738717E-4+t33*t57*t72*2.848666080295449E-2-t33*t58*t71*4.353713116283091E-4-t56*t57*t58*5.585324801683264E-2-t33*t71*t72*4.353713116283091E-4-t56*t57*t72*5.585324801683264E-2+t56*t58*t71*5.568788485803049E-2+t57*t58*t70*5.585324801683264E-2+t56*t71*t72*5.568788485803049E-2+t57*t70*t72*5.585324801683264E-2-t58*t70*t71*5.568788485803049E-2-t70*t71*t72*5.568788485803049E-2-t9*t10*t33*t34*2.230820224171545E-1+t9*t10*t34*t58*2.896401859571867E-2+t9*t10*t34*t72*2.896401859571867E-2+t10*t33*t34*t58*1.113095897646181E-3+t10*t33*t34*t72*1.113095897646181E-3+t9*t33*t57*t58*3.159934013264157E-2+t10*t34*t56*t58*3.156734948660587E-2+t9*t33*t57*t72*3.159934013264157E-2-t9*t33*t58*t71*3.153538479057603E-2+t10*t34*t56*t72*3.156734948660587E-2-t10*t34*t58*t70*3.156734948660587E-2-t9*t33*t71*t72*3.153538479057603E-2-t10*t34*t70*t72*3.156734948660587E-2+t9*t10*t33*t34*t58*2.230820224171545E-1+t9*t10*t33*t34*t72*2.230820224171545E-1;
    J(10, 17) = t57*4.176918863775431E-6-t71*4.176918863775431E-6+t10*t34*9.830160358935766E-5+t9*t57*1.113095897646181E-3-t9*t71*1.113095897646181E-3-t33*t57*2.896401859571867E-2+t33*t71*2.896401859571867E-2+t56*t57*5.576841714485674E-2-t57*t58*4.176918863775431E-6-t56*t71*5.576841714485674E-2-t57*t70*5.577259406372051E-2-t57*t72*4.176918863775431E-6+t58*t71*4.176918863775431E-6+t70*t71*5.577259406372051E-2+t71*t72*4.176918863775431E-6+t9*t10*t11*2.87566285868891E-1-t9*t10*t34*5.784406422916559E-2+t9*t11*t34*3.104080344633556E-3+t10*t11*t33*3.064210757541298E-3-t10*t33*t34*2.220487864525553E-3-t11*t33*t34*2.879825211612492E-1-t9*t33*t57*3.156734948660587E-2-t10*t34*t56*3.161651326340348E-2-t10*t34*t58*9.830160358935766E-5+t9*t33*t71*3.156734948660587E-2+t10*t34*t70*3.151821165981412E-2-t10*t34*t72*9.830160358935766E-5-t9*t57*t58*1.113095897646181E-3-t9*t57*t72*1.113095897646181E-3+t9*t58*t71*1.113095897646181E-3+t33*t57*t58*2.896401859571867E-2+t9*t71*t72*1.113095897646181E-3+t33*t57*t72*2.896401859571867E-2-t33*t58*t71*2.896401859571867E-2-t56*t57*t58*5.576841714485674E-2-t33*t71*t72*2.896401859571867E-2-t56*t57*t72*5.576841714485674E-2+t56*t58*t71*5.576841714485674E-2+t57*t58*t70*5.577259406372051E-2+t56*t71*t72*5.576841714485674E-2+t57*t70*t72*5.577259406372051E-2-t58*t70*t71*5.577259406372051E-2-t70*t71*t72*5.577259406372051E-2-t9*t10*t33*t34*2.230822657497263E-1+t9*t10*t34*t58*5.784406422916559E-2+t9*t10*t34*t72*5.784406422916559E-2+t10*t33*t34*t58*2.220487864525553E-3+t10*t33*t34*t72*2.220487864525553E-3+t9*t33*t57*t58*3.156734948660587E-2+t10*t34*t56*t58*3.161651326340348E-2+t9*t33*t57*t72*3.156734948660587E-2-t9*t33*t58*t71*3.156734948660587E-2+t10*t34*t56*t72*3.161651326340348E-2-t10*t34*t58*t70*3.151821165981412E-2-t9*t33*t71*t72*3.156734948660587E-2-t10*t34*t70*t72*3.151821165981412E-2+t9*t10*t33*t34*t58*2.230822657497263E-1+t9*t10*t33*t34*t72*2.230822657497263E-1;
    J(10, 18) = t9*t10*t35*3.104080344633556E-3-t9*t11*t57*1.549180159908108E-2-t9*t34*t35*2.87566285868891E-1-t10*t33*t35*2.879825211612492E-1-t9*t11*t71*1.549180159908108E-2-t11*t33*t57*8.234792310892687E-4-t33*t34*t35*3.064210757541298E-3-t11*t33*t71*8.234792310892687E-4;
    J(10, 19) = t85+t104+t118+t120+t122+t132+t134+t136+t137+t138+t141+t152+t153-t12*t13*9.070578524442012E-3-t12*t37*1.718955829676478E-2-t13*t36*2.991020934719677E-3-t36*t37*5.668252425759201E-3;
    J(10, 20) = t12*t13*5.668252425759201E-3-t12*t37*2.991020934719677E-3-t13*t36*1.718955829676478E-2+t36*t37*9.070578524442012E-3;
    J(11, 16) = t9*t57*-2.99216464372371E-3+t9*t58*3.029602411709927E-3-t9*t71*3.743776798621798E-5+t9*t72*3.029602411709927E-3+t33*t57*8.817574215578565E-4-t33*t58*7.546580486902988E-4-t33*t71*1.270993728675576E-4-t33*t72*7.546580486902988E-4+t9*t10*t34*1.009577116650697E-3+t10*t33*t34*2.959208971395655E-3-t9*t57*t58*3.743776798621798E-5-t9*t57*t72*3.743776798621798E-5-t9*t58*t71*2.99216464372371E-3-t33*t57*t58*1.270993728675576E-4-t9*t71*t72*2.99216464372371E-3-t33*t57*t72*1.270993728675576E-4+t33*t58*t71*8.817574215578565E-4+t33*t71*t72*8.817574215578565E-4-t9*t10*t34*t58*1.009577116650697E-3-t9*t10*t34*t72*1.009577116650697E-3-t10*t33*t34*t58*2.959208971395655E-3-t10*t33*t34*t72*2.959208971395655E-3;
    J(11, 17) = t9*t57*-2.959208971395655E-3+t9*t71*2.959208971395655E-3+t33*t57*1.009577116650697E-3-t33*t71*1.009577116650697E-3+t56*t57*2.297574320158763E-5-t56*t71*2.297574320158763E-5+t57*t70*2.297574320158763E-5-t70*t71*2.297574320158763E-5+t9*t10*t34*2.017713588850828E-3-t10*t11*t56*1.549973690114125E-2+t10*t33*t34*5.909453751474983E-3-t10*t11*t70*1.549973690114125E-2+t10*t34*t56*3.122371117714521E-3+t11*t34*t56*6.574122182670539E-4+t10*t34*t70*3.122371117714521E-3+t11*t34*t70*6.574122182670539E-4+t9*t57*t58*2.959208971395655E-3+t9*t57*t72*2.959208971395655E-3-t9*t58*t71*2.959208971395655E-3-t33*t57*t58*1.009577116650697E-3-t9*t71*t72*2.959208971395655E-3-t33*t57*t72*1.009577116650697E-3+t33*t58*t71*1.009577116650697E-3-t56*t57*t58*2.297574320158763E-5+t33*t71*t72*1.009577116650697E-3-t56*t57*t72*2.297574320158763E-5+t56*t58*t71*2.297574320158763E-5-t57*t58*t70*2.297574320158763E-5+t56*t71*t72*2.297574320158763E-5-t57*t70*t72*2.297574320158763E-5+t58*t70*t71*2.297574320158763E-5+t70*t71*t72*2.297574320158763E-5-t9*t10*t34*t58*2.017713588850828E-3-t9*t10*t34*t72*2.017713588850828E-3-t10*t33*t34*t58*5.909453751474983E-3-t10*t33*t34*t72*5.909453751474983E-3-t10*t34*t56*t58*3.122371117714521E-3-t10*t34*t56*t72*3.122371117714521E-3-t10*t34*t58*t70*3.122371117714521E-3-t10*t34*t70*t72*3.122371117714521E-3;
    J(11, 18) = t10*t35*t56*6.574122182670539E-4+t10*t35*t70*6.574122182670539E-4-t11*t56*t57*2.875818595898751E-1+t34*t35*t56*1.549973690114125E-2-t11*t56*t71*2.875818595898751E-1-t11*t57*t70*2.875818595898751E-1+t34*t35*t70*1.549973690114125E-2-t11*t70*t71*2.875818595898751E-1;
    J(11, 20) = t121+t125+t37*t59*1.81E-2+t37*t73*1.81E-2;
    J(12, 24) = t61*1.586414893789752E-3-t75*1.586414893789752E-3-t14*t15*4.954563135453205E-1+t14*t38*8.359461048286212E-4-t15*t38*6.725214316796237E-2-t61*t62*1.586414893789752E-3-t61*t76*1.586414893789752E-3+t62*t75*1.586414893789752E-3+t75*t76*1.586414893789752E-3-t14*t38*t62*8.359461048286212E-4-t14*t38*t76*8.359461048286212E-4;
    J(12, 25) = t14*t39*-6.725214316796237E-2-t15*t61*4.365860704565494E-4+t38*t39*4.954563135453205E-1-t15*t75*4.365860704565494E-4;
    J(12, 26) = t16*2.47718E-1+t40*3.370300000000001E-2+t90+t100+t110+t112+t16*t63*2.47718E-1+t16*t77*2.47718E-1+t40*t63*3.370300000000001E-2+t40*t77*3.3703E-2;
    J(12, 27) = t90+t100+t110+t112+t16*t17*t41*1.387778780781446E-17-t17*t40*t41*5.551115123125783E-17;
    J(13, 24) = t61*-4.170075425381788E-4+t75*4.170075425381788E-4+t14*t15*6.711175107535902E-2+t14*t38*3.164698577274911E-3-t15*t38*4.941905177865888E-1+t61*t62*4.170075425381788E-4+t61*t76*4.170075425381788E-4-t62*t75*4.170075425381788E-4-t75*t76*4.170075425381788E-4-t14*t38*t62*3.164698577274911E-3-t14*t38*t76*3.164698577274911E-3;
    J(13, 25) = t14*t39*-4.941905177865888E-1-t15*t61*3.566153386244494E-2-t38*t39*6.711175107535902E-2-t15*t75*3.566153386244494E-2;
    J(13, 26) = t16*-3.370300000000001E-2+t40*2.47718E-1+t94+t96+t103+t116-t16*t63*3.370300000000001E-2-t16*t77*3.3703E-2+t40*t63*2.47718E-1+t40*t77*2.47718E-1;
    J(13, 27) = t94+t96+t103+t116+t16*t17*t41*5.551115123125783E-17+t17*t40*t41*1.387778780781446E-17;
    J(14, 24) = t61*2.84294570084937E-5-t75*2.84294570084937E-5-t14*t15*4.365115769377904E-3-t14*t38*2.270241925564839E-4+t15*t38*3.539606431708408E-2-t61*t62*2.84294570084937E-5-t61*t76*2.84294570084937E-5+t62*t75*2.84294570084937E-5+t75*t76*2.84294570084937E-5+t14*t38*t62*2.270241925564839E-4+t14*t38*t76*2.270241925564839E-4;
    J(14, 25) = t14*t39*3.539606431708408E-2-t15*t61*4.987264424463382E-1+t38*t39*4.365115769377904E-3-t15*t75*4.987264424463382E-1;
    J(15, 9) = t50*-1.586414893789752E-3+t64*1.586414893789752E-3+t2*t3*4.954563135453205E-1+t2*t26*8.359461048286212E-4-t3*t26*6.725214316796237E-2+t50*t51*1.586414893789752E-3+t50*t65*1.586414893789752E-3-t51*t64*1.586414893789752E-3-t64*t65*1.586414893789752E-3-t2*t26*t51*8.359461048286212E-4-t2*t26*t65*8.359461048286212E-4;
    J(15, 10) = t2*t27*-6.725214316796237E-2-t3*t50*4.365860704565494E-4-t26*t27*4.954563135453205E-1-t3*t64*4.365860704565494E-4;
    J(15, 11) = t4*-2.47718E-1+t28*3.370300000000001E-2+t89+t98+t106+t108-t4*t52*2.47718E-1-t4*t66*2.47718E-1+t28*t52*3.370300000000001E-2+t28*t66*3.3703E-2;
    J(15, 12) = t89+t98+t106+t108+t4*t5*t29*1.387778780781446E-17+t5*t28*t29*5.551115123125783E-17;
    J(16, 9) = t50*-4.170075425381788E-4+t64*4.170075425381788E-4+t2*t3*6.711175107535902E-2-t2*t26*3.164698577274911E-3+t3*t26*4.941905177865888E-1+t50*t51*4.170075425381788E-4+t50*t65*4.170075425381788E-4-t51*t64*4.170075425381788E-4-t64*t65*4.170075425381788E-4+t2*t26*t51*3.164698577274911E-3+t2*t26*t65*3.164698577274911E-3;
    J(16, 10) = t2*t27*4.941905177865888E-1+t3*t50*3.566153386244494E-2-t26*t27*6.711175107535902E-2+t3*t64*3.566153386244494E-2;
    J(16, 11) = t4*-3.370300000000001E-2-t28*2.47718E-1+t92+t93+t102+t114-t4*t52*3.370300000000001E-2-t4*t66*3.3703E-2-t28*t52*2.47718E-1-t28*t66*2.47718E-1;
    J(16, 12) = t92+t93+t102+t114-t4*t5*t29*5.551115123125783E-17+t5*t28*t29*1.387778780781446E-17;
    J(17, 9) = t50*-2.84294570084937E-5+t64*2.84294570084937E-5+t2*t3*4.365115769377904E-3-t2*t26*2.270241925564839E-4+t3*t26*3.539606431708408E-2+t50*t51*2.84294570084937E-5+t50*t65*2.84294570084937E-5-t51*t64*2.84294570084937E-5-t64*t65*2.84294570084937E-5+t2*t26*t51*2.270241925564839E-4+t2*t26*t65*2.270241925564839E-4;
    J(17, 10) = t2*t27*3.539606431708408E-2-t3*t50*4.987264424463382E-1-t26*t27*4.365115769377904E-3-t3*t64*4.987264424463382E-1;

    fkPtr_->compute(0, 
                    contact_joint_id, 
                    q, 
                    nullptr, 
                    &stance_foot_endT, 
                    1);
                    
    J.block(18, 0, 3, modelPtr_->nv) = fkPtr_->getTranslationJacobian();
    J.block(21, 0, 3, modelPtr_->nv) = fkPtr_->getRPYJacobian();
}

void DigitDynamicsConstraints::get_Jx_partial_dq(const VecX& q, const VecX& x) {
    assert(x.size() == modelPtr_->nv);
    assert(Jx_partial_dq.rows() == NUM_DEPENDENT_JOINTS);
    assert(Jx_partial_dq.cols() == modelPtr_->nv);

    double x1 = x(0);
    double x2 = x(1);
    double x3 = x(2);
    double x4 = x(3);
    double x5 = x(4);
    double x6 = x(5);
    double x7 = x(6);
    double x8 = x(7);
    double x9 = x(8);
    double x10 = x(9);
    double x11 = x(10);
    double x12 = x(11);
    double x13 = x(12);
    double x14 = x(13);
    double x15 = x(14);
    double x16 = x(15);
    double x17 = x(16);
    double x18 = x(17);
    double x19 = x(18);
    double x20 = x(19);
    double x21 = x(20);
    double x22 = x(21);
    double x23 = x(22);
    double x24 = x(23);
    double x25 = x(24);
    double x26 = x(25);
    double x27 = x(26);
    double x28 = x(27);
    double x29 = x(28);
    double x30 = x(29);
    double x31 = x(30);
    double x32 = x(31);
    double x33 = x(32);
    double x34 = x(33);
    double x35 = x(34);
    double x36 = x(35);

    double t2 = cos(q(9));
    double t3 = cos(q(10));
    double t4 = cos(q(11));
    double t5 = cos(q(12));
    double t6 = cos(q(13));
    double t7 = cos(q(14));
    double t8 = cos(q(15));
    double t9 = cos(q(16));
    double t10 = cos(q(17));
    double t11 = cos(q(18));
    double t12 = cos(q(19));
    double t13 = cos(q(20));
    double t14 = cos(q(24));
    double t15 = cos(q(25));
    double t16 = cos(q(26));
    double t17 = cos(q(27));
    double t18 = cos(q(28));
    double t19 = cos(q(29));
    double t20 = cos(q(30));
    double t21 = cos(q(31));
    double t22 = cos(q(32));
    double t23 = cos(q(33));
    double t24 = cos(q(34));
    double t25 = cos(q(35));
    double t26 = sin(q(9));
    double t27 = sin(q(10));
    double t28 = sin(q(11));
    double t29 = sin(q(12));
    double t30 = sin(q(13));
    double t31 = sin(q(14));
    double t32 = sin(q(15));
    double t33 = sin(q(16));
    double t34 = sin(q(17));
    double t35 = sin(q(18));
    double t36 = sin(q(19));
    double t37 = sin(q(20));
    double t38 = sin(q(24));
    double t39 = sin(q(25));
    double t40 = sin(q(26));
    double t41 = sin(q(27));
    double t42 = sin(q(28));
    double t43 = sin(q(29));
    double t44 = sin(q(30));
    double t45 = sin(q(31));
    double t46 = sin(q(32));
    double t47 = sin(q(33));
    double t48 = sin(q(34));
    double t49 = sin(q(35));
    double t50 = t12*1.696216709330505E-2;
    double t51 = t24*1.696216709330505E-2;
    double t52 = t36*1.696216709330505E-2;
    double t53 = t48*-1.696216709330505E-2;
    double t54 = t48*1.696216709330505E-2;
    double t55 = t37*9.551E-3;
    double t56 = t49*9.551E-3;
    double t57 = t12*5.143951577823025E-2;
    double t58 = t24*5.143951577823025E-2;
    double t59 = t36*5.143951577823025E-2;
    double t60 = t48*5.143951577823025E-2;
    double t61 = t4*t5*-6.370345907739961E-5;
    double t62 = t4*t5*6.370345907739961E-5;
    double t63 = t16*t17*6.370345907739961E-5;
    double t64 = t2*t3*6.725214316796237E-2;
    double t65 = t14*t15*6.725214316796237E-2;
    double t66 = t12*t13*1.69996184260823E-2;
    double t67 = t24*t25*1.69996184260823E-2;
    double t68 = t12*t13*1.718955829676478E-2;
    double t69 = t24*t25*1.718955829676478E-2;
    double t70 = t2*t3*3.539606431708408E-2;
    double t71 = t14*t15*3.539606431708408E-2;
    double t72 = t6*t8*-2.246221860400801E-3;
    double t73 = t6*t8*2.246221860400801E-3;
    double t74 = t18*t20*2.246221860400801E-3;
    double t75 = t12*t13*9.070578524442012E-3;
    double t76 = t12*t13*9.070578524442013E-3;
    double t77 = t24*t25*9.070578524442012E-3;
    double t78 = t24*t25*9.070578524442013E-3;
    double t79 = t7*t8*2.360206106172353E-3;
    double t80 = t19*t20*2.360206106172353E-3;
    double t81 = t4*t29*6.370345907739961E-5;
    double t82 = t5*t28*6.370345907739961E-5;
    double t83 = t16*t41*6.370345907739961E-5;
    double t84 = t17*t40*6.370345907739961E-5;
    double t85 = t2*t27*6.711175107535902E-2;
    double t86 = t3*t26*6.711175107535902E-2;
    double t87 = t14*t39*6.711175107535902E-2;
    double t88 = t15*t38*6.711175107535902E-2;
    double t89 = t12*t37*1.69996184260823E-2;
    double t90 = t13*t36*-1.69996184260823E-2;
    double t91 = t13*t36*1.69996184260823E-2;
    double t92 = t24*t49*1.69996184260823E-2;
    double t93 = t25*t48*1.69996184260823E-2;
    double t94 = t10*t11*6.574122182670539E-4;
    double t95 = t22*t23*6.574122182670539E-4;
    double t96 = t12*t37*1.718955829676478E-2;
    double t97 = t13*t36*-1.718955829676478E-2;
    double t98 = t13*t36*1.718955829676478E-2;
    double t99 = t24*t49*1.718955829676478E-2;
    double t100 = t25*t48*1.718955829676478E-2;
    double t101 = t2*t27*4.365115769377904E-3;
    double t102 = t3*t26*-4.365115769377904E-3;
    double t103 = t3*t26*4.365115769377904E-3;
    double t104 = t14*t39*4.365115769377904E-3;
    double t105 = t15*t38*4.365115769377904E-3;
    double t106 = t6*t32*2.246221860400801E-3;
    double t107 = t8*t30*2.246221860400801E-3;
    double t108 = t18*t44*2.246221860400801E-3;
    double t109 = t20*t42*2.246221860400801E-3;
    double t110 = t12*t37*-9.070578524442012E-3;
    double t111 = t12*t37*9.070578524442012E-3;
    double t112 = t12*t37*9.070578524442013E-3;
    double t113 = t13*t36*-9.070578524442013E-3;
    double t114 = t13*t36*-9.070578524442012E-3;
    double t115 = t13*t36*9.070578524442012E-3;
    double t116 = t13*t36*9.070578524442013E-3;
    double t117 = t24*t49*-9.070578524442013E-3;
    double t118 = t24*t49*9.070578524442012E-3;
    double t119 = t24*t49*9.070578524442013E-3;
    double t120 = t25*t48*9.070578524442012E-3;
    double t121 = t25*t48*9.070578524442013E-3;
    double t122 = t12*t13*5.605619802270151E-3;
    double t123 = t24*t25*5.605619802270151E-3;
    double t124 = t12*t13*5.668252425759201E-3;
    double t125 = t24*t25*5.668252425759201E-3;
    double t126 = t12*t13*2.991020934719675E-3;
    double t127 = t12*t13*2.991020934719677E-3;
    double t128 = t24*t25*2.991020934719675E-3;
    double t129 = t24*t25*2.991020934719677E-3;
    double t130 = t28*t29*6.370345907739961E-5;
    double t131 = t40*t41*-6.370345907739961E-5;
    double t132 = t40*t41*6.370345907739961E-5;
    double t133 = t26*t27*-6.725214316796237E-2;
    double t134 = t26*t27*6.725214316796237E-2;
    double t135 = t38*t39*6.725214316796237E-2;
    double t136 = t36*t37*1.69996184260823E-2;
    double t137 = t48*t49*1.69996184260823E-2;
    double t138 = t36*t37*-1.718955829676478E-2;
    double t139 = t36*t37*1.718955829676478E-2;
    double t140 = t48*t49*-1.718955829676478E-2;
    double t141 = t48*t49*1.718955829676478E-2;
    double t142 = t6*t8*1.26575449899492E-2;
    double t143 = t18*t20*-1.26575449899492E-2;
    double t144 = t18*t20*1.26575449899492E-2;
    double t145 = t26*t27*3.539606431708408E-2;
    double t146 = t38*t39*-3.539606431708408E-2;
    double t147 = t38*t39*3.539606431708408E-2;
    double t148 = t30*t32*2.246221860400801E-3;
    double t149 = t42*t44*-2.246221860400801E-3;
    double t150 = t42*t44*2.246221860400801E-3;
    double t151 = t36*t37*-9.070578524442013E-3;
    double t152 = t36*t37*9.070578524442012E-3;
    double t153 = t36*t37*9.070578524442013E-3;
    double t154 = t48*t49*-9.070578524442013E-3;
    double t155 = t48*t49*9.070578524442012E-3;
    double t156 = t48*t49*9.070578524442013E-3;
    double t157 = t12*t37*-5.605619802270151E-3;
    double t158 = t12*t37*5.605619802270151E-3;
    double t159 = t13*t36*5.605619802270151E-3;
    double t160 = t24*t49*5.605619802270151E-3;
    double t161 = t25*t48*-5.605619802270151E-3;
    double t162 = t25*t48*5.605619802270151E-3;
    double t163 = t9*t11*8.234792310892687E-4;
    double t164 = t21*t23*8.234792310892687E-4;
    double t165 = t12*t37*5.668252425759201E-3;
    double t166 = t13*t36*5.668252425759201E-3;
    double t167 = t24*t49*-5.668252425759201E-3;
    double t168 = t24*t49*5.668252425759201E-3;
    double t169 = t25*t48*-5.668252425759201E-3;
    double t170 = t25*t48*5.668252425759201E-3;
    double t171 = t31*t32*-2.360206106172353E-3;
    double t172 = t31*t32*2.360206106172353E-3;
    double t173 = t43*t44*2.360206106172353E-3;
    double t174 = t12*t37*-2.991020934719677E-3;
    double t175 = t12*t37*2.991020934719675E-3;
    double t176 = t12*t37*2.991020934719677E-3;
    double t177 = t13*t36*2.991020934719675E-3;
    double t178 = t13*t36*2.991020934719677E-3;
    double t179 = t24*t49*-2.991020934719675E-3;
    double t180 = t24*t49*2.991020934719675E-3;
    double t181 = t24*t49*2.991020934719677E-3;
    double t182 = t25*t48*2.991020934719675E-3;
    double t183 = t25*t48*2.991020934719677E-3;
    double t184 = t34*t35*-6.574122182670539E-4;
    double t185 = t34*t35*6.574122182670539E-4;
    double t186 = t46*t47*6.574122182670539E-4;
    double t187 = t7*t32*1.263678697118524E-2;
    double t188 = t8*t31*1.263678697118524E-2;
    double t189 = t19*t44*1.263678697118524E-2;
    double t190 = t20*t43*-1.263678697118524E-2;
    double t191 = t20*t43*1.263678697118524E-2;
    double t192 = t6*t32*1.26575449899492E-2;
    double t193 = t8*t30*1.26575449899492E-2;
    double t194 = t18*t44*1.26575449899492E-2;
    double t195 = t20*t42*1.26575449899492E-2;
    double t196 = t4*t5*1.20005624647208E-1;
    double t197 = t16*t17*-1.20005624647208E-1;
    double t198 = t16*t17*1.20005624647208E-1;
    double t199 = t36*t37*5.605619802270151E-3;
    double t200 = t48*t49*5.605619802270151E-3;
    double t201 = t9*t35*-8.234792310892687E-4;
    double t202 = t9*t35*8.234792310892687E-4;
    double t203 = t11*t33*8.234792310892687E-4;
    double t204 = t21*t47*8.234792310892687E-4;
    double t205 = t23*t45*8.234792310892687E-4;
    double t206 = t36*t37*5.668252425759201E-3;
    double t207 = t48*t49*5.668252425759201E-3;
    double t208 = t2*t3*4.941905177865888E-1;
    double t209 = t14*t15*4.941905177865888E-1;
    double t210 = t9*t11*1.549180159908108E-2;
    double t211 = t21*t23*-1.549180159908108E-2;
    double t212 = t21*t23*1.549180159908108E-2;
    double t213 = t36*t37*-2.991020934719677E-3;
    double t214 = t36*t37*2.991020934719675E-3;
    double t215 = t36*t37*2.991020934719677E-3;
    double t216 = t48*t49*-2.991020934719677E-3;
    double t217 = t48*t49*2.991020934719675E-3;
    double t218 = t48*t49*2.991020934719677E-3;
    double t219 = t30*t32*-1.26575449899492E-2;
    double t220 = t30*t32*1.26575449899492E-2;
    double t221 = t42*t44*1.26575449899492E-2;
    double t222 = t4*t29*1.20005624647208E-1;
    double t223 = t5*t28*1.20005624647208E-1;
    double t224 = t16*t41*1.20005624647208E-1;
    double t225 = t17*t40*1.20005624647208E-1;
    double t226 = t33*t35*8.234792310892687E-4;
    double t227 = t45*t47*-8.234792310892687E-4;
    double t228 = t45*t47*8.234792310892687E-4;
    double t229 = t2*t27*4.954563135453205E-1;
    double t230 = t3*t26*4.954563135453205E-1;
    double t231 = t14*t39*4.954563135453205E-1;
    double t232 = t15*t38*4.954563135453205E-1;
    double t233 = t9*t35*1.549180159908108E-2;
    double t234 = t11*t33*-1.549180159908108E-2;
    double t235 = t11*t33*1.549180159908108E-2;
    double t236 = t21*t47*1.549180159908108E-2;
    double t237 = t23*t45*1.549180159908108E-2;
    double t238 = t10*t35*1.549973690114125E-2;
    double t239 = t11*t34*1.549973690114125E-2;
    double t240 = t22*t47*1.549973690114125E-2;
    double t241 = t23*t46*-1.549973690114125E-2;
    double t242 = t23*t46*1.549973690114125E-2;
    double t243 = t28*t29*-1.20005624647208E-1;
    double t244 = t28*t29*1.20005624647208E-1;
    double t245 = t40*t41*1.20005624647208E-1;
    double t246 = t26*t27*4.941905177865888E-1;
    double t247 = t38*t39*-4.941905177865888E-1;
    double t248 = t38*t39*4.941905177865888E-1;
    double t249 = t33*t35*1.549180159908108E-2;
    double t250 = t45*t47*1.549180159908108E-2;
    double t251 = t9*t10*t11*2.87566285868891E-1;
    double t252 = t21*t22*t23*-2.87566285868891E-1;
    double t253 = t21*t22*t23*2.87566285868891E-1;
    double t254 = t9*t10*t11*2.879825211612492E-1;
    double t255 = t21*t22*t23*-2.879825211612492E-1;
    double t256 = t21*t22*t23*2.879825211612492E-1;
    double t257 = t6*t7*t8*-3.397508858570615E-1;
    double t258 = t6*t7*t8*3.397508858570615E-1;
    double t259 = t18*t19*t20*3.397508858570615E-1;
    double t260 = t6*t7*t8*-3.399783924207052E-1;
    double t261 = t6*t7*t8*3.399783924207052E-1;
    double t262 = t18*t19*t20*3.399783924207052E-1;
    double t263 = t9*t10*t35*2.87566285868891E-1;
    double t264 = t9*t11*t34*-2.87566285868891E-1;
    double t265 = t9*t11*t34*2.87566285868891E-1;
    double t266 = t10*t11*t33*2.87566285868891E-1;
    double t267 = t21*t22*t47*-2.87566285868891E-1;
    double t268 = t21*t22*t47*2.87566285868891E-1;
    double t269 = t21*t23*t46*2.87566285868891E-1;
    double t270 = t22*t23*t45*2.87566285868891E-1;
    double t271 = t9*t10*t35*2.879825211612492E-1;
    double t272 = t9*t11*t34*2.879825211612492E-1;
    double t273 = t10*t11*t33*-2.879825211612492E-1;
    double t274 = t10*t11*t33*2.879825211612492E-1;
    double t275 = t21*t22*t47*2.879825211612492E-1;
    double t276 = t21*t23*t46*2.879825211612492E-1;
    double t277 = t22*t23*t45*2.879825211612492E-1;
    double t278 = t6*t7*t8*3.020283789547073E-3;
    double t279 = t18*t19*t20*-3.020283789547073E-3;
    double t280 = t18*t19*t20*3.020283789547073E-3;
    double t281 = t9*t10*t11*-3.064210757541298E-3;
    double t282 = t9*t10*t11*3.064210757541298E-3;
    double t283 = t21*t22*t23*3.064210757541298E-3;
    double t284 = t9*t10*t11*3.104080344633556E-3;
    double t285 = t21*t22*t23*3.104080344633556E-3;
    double t286 = t6*t7*t8*3.105990081579729E-3;
    double t287 = t18*t19*t20*-3.105990081579729E-3;
    double t288 = t18*t19*t20*3.105990081579729E-3;
    double t289 = t6*t7*t32*-3.397508858570615E-1;
    double t290 = t6*t7*t32*3.397508858570615E-1;
    double t291 = t6*t8*t31*3.397508858570615E-1;
    double t292 = t7*t8*t30*3.397508858570615E-1;
    double t293 = t18*t19*t44*3.397508858570615E-1;
    double t294 = t18*t20*t43*3.397508858570615E-1;
    double t295 = t19*t20*t42*3.397508858570615E-1;
    double t296 = t6*t7*t32*-3.399783924207052E-1;
    double t297 = t6*t7*t32*3.399783924207052E-1;
    double t298 = t6*t8*t31*3.399783924207052E-1;
    double t299 = t7*t8*t30*3.399783924207052E-1;
    double t300 = t18*t19*t44*3.399783924207052E-1;
    double t301 = t18*t20*t43*3.399783924207052E-1;
    double t302 = t19*t20*t42*3.399783924207052E-1;
    double t303 = t9*t34*t35*2.87566285868891E-1;
    double t304 = t10*t33*t35*2.87566285868891E-1;
    double t305 = t11*t33*t34*-2.87566285868891E-1;
    double t306 = t11*t33*t34*2.87566285868891E-1;
    double t307 = t21*t46*t47*2.87566285868891E-1;
    double t308 = t22*t45*t47*2.87566285868891E-1;
    double t309 = t23*t45*t46*2.87566285868891E-1;
    double t310 = t9*t34*t35*2.879825211612492E-1;
    double t311 = t10*t33*t35*2.879825211612492E-1;
    double t312 = t11*t33*t34*-2.879825211612492E-1;
    double t313 = t11*t33*t34*2.879825211612492E-1;
    double t314 = t21*t46*t47*2.879825211612492E-1;
    double t315 = t22*t45*t47*2.879825211612492E-1;
    double t316 = t23*t45*t46*2.879825211612492E-1;
    double t317 = t6*t7*t32*3.020283789547073E-3;
    double t318 = t6*t8*t31*3.020283789547073E-3;
    double t319 = t7*t8*t30*3.020283789547073E-3;
    double t320 = t18*t19*t44*-3.020283789547073E-3;
    double t321 = t18*t19*t44*3.020283789547073E-3;
    double t322 = t18*t20*t43*3.020283789547073E-3;
    double t323 = t19*t20*t42*3.020283789547073E-3;
    double t324 = t9*t10*t35*-3.064210757541298E-3;
    double t325 = t9*t10*t35*3.064210757541298E-3;
    double t326 = t9*t11*t34*3.064210757541298E-3;
    double t327 = t10*t11*t33*3.064210757541298E-3;
    double t328 = t21*t22*t47*3.064210757541298E-3;
    double t329 = t21*t23*t46*3.064210757541298E-3;
    double t330 = t22*t23*t45*3.064210757541298E-3;
    double t331 = t9*t10*t35*-3.104080344633556E-3;
    double t332 = t9*t10*t35*3.104080344633556E-3;
    double t333 = t9*t11*t34*3.104080344633556E-3;
    double t334 = t10*t11*t33*3.104080344633556E-3;
    double t335 = t21*t22*t47*3.104080344633556E-3;
    double t336 = t21*t23*t46*3.104080344633556E-3;
    double t337 = t22*t23*t45*3.104080344633556E-3;
    double t338 = t6*t7*t32*3.105990081579729E-3;
    double t339 = t6*t8*t31*3.105990081579729E-3;
    double t340 = t7*t8*t30*3.105990081579729E-3;
    double t341 = t18*t19*t44*-3.105990081579729E-3;
    double t342 = t18*t19*t44*3.105990081579729E-3;
    double t343 = t18*t20*t43*3.105990081579729E-3;
    double t344 = t19*t20*t42*3.105990081579729E-3;
    double t345 = t6*t31*t32*3.397508858570615E-1;
    double t346 = t7*t30*t32*3.397508858570615E-1;
    double t347 = t8*t30*t31*3.397508858570615E-1;
    double t348 = t18*t43*t44*3.397508858570615E-1;
    double t349 = t19*t42*t44*3.397508858570615E-1;
    double t350 = t20*t42*t43*-3.397508858570615E-1;
    double t351 = t20*t42*t43*3.397508858570615E-1;
    double t352 = t6*t31*t32*3.399783924207052E-1;
    double t353 = t7*t30*t32*3.399783924207052E-1;
    double t354 = t8*t30*t31*3.399783924207052E-1;
    double t355 = t18*t43*t44*3.399783924207052E-1;
    double t356 = t19*t42*t44*3.399783924207052E-1;
    double t357 = t20*t42*t43*-3.399783924207052E-1;
    double t358 = t20*t42*t43*3.399783924207052E-1;
    double t359 = t33*t34*t35*-2.87566285868891E-1;
    double t360 = t33*t34*t35*2.87566285868891E-1;
    double t361 = t45*t46*t47*-2.87566285868891E-1;
    double t362 = t45*t46*t47*2.87566285868891E-1;
    double t363 = t33*t34*t35*-2.879825211612492E-1;
    double t364 = t33*t34*t35*2.879825211612492E-1;
    double t365 = t45*t46*t47*2.879825211612492E-1;
    double t366 = t6*t31*t32*3.020283789547073E-3;
    double t367 = t7*t30*t32*3.020283789547073E-3;
    double t368 = t8*t30*t31*-3.020283789547073E-3;
    double t369 = t8*t30*t31*3.020283789547073E-3;
    double t370 = t18*t43*t44*3.020283789547073E-3;
    double t371 = t19*t42*t44*3.020283789547073E-3;
    double t372 = t20*t42*t43*3.020283789547073E-3;
    double t373 = t9*t34*t35*3.064210757541298E-3;
    double t374 = t10*t33*t35*3.064210757541298E-3;
    double t375 = t11*t33*t34*-3.064210757541298E-3;
    double t376 = t11*t33*t34*3.064210757541298E-3;
    double t377 = t21*t46*t47*-3.064210757541298E-3;
    double t378 = t21*t46*t47*3.064210757541298E-3;
    double t379 = t22*t45*t47*3.064210757541298E-3;
    double t380 = t23*t45*t46*-3.064210757541298E-3;
    double t381 = t23*t45*t46*3.064210757541298E-3;
    double t382 = t9*t34*t35*3.104080344633556E-3;
    double t383 = t10*t33*t35*3.104080344633556E-3;
    double t384 = t11*t33*t34*3.104080344633556E-3;
    double t385 = t21*t46*t47*3.104080344633556E-3;
    double t386 = t22*t45*t47*-3.104080344633556E-3;
    double t387 = t22*t45*t47*3.104080344633556E-3;
    double t388 = t23*t45*t46*-3.104080344633556E-3;
    double t389 = t23*t45*t46*3.104080344633556E-3;
    double t390 = t6*t31*t32*3.105990081579729E-3;
    double t391 = t7*t30*t32*3.105990081579729E-3;
    double t392 = t8*t30*t31*-3.105990081579729E-3;
    double t393 = t8*t30*t31*3.105990081579729E-3;
    double t394 = t18*t43*t44*3.105990081579729E-3;
    double t395 = t19*t42*t44*3.105990081579729E-3;
    double t396 = t20*t42*t43*3.105990081579729E-3;
    double t397 = t30*t31*t32*3.397508858570615E-1;
    double t398 = t42*t43*t44*-3.397508858570615E-1;
    double t399 = t42*t43*t44*3.397508858570615E-1;
    double t400 = t30*t31*t32*3.399783924207052E-1;
    double t401 = t42*t43*t44*-3.399783924207052E-1;
    double t402 = t42*t43*t44*3.399783924207052E-1;
    double t403 = t30*t31*t32*-3.020283789547073E-3;
    double t404 = t30*t31*t32*3.020283789547073E-3;
    double t405 = t42*t43*t44*3.020283789547073E-3;
    double t406 = t33*t34*t35*3.064210757541298E-3;
    double t407 = t45*t46*t47*-3.064210757541298E-3;
    double t408 = t45*t46*t47*3.064210757541298E-3;
    double t409 = t33*t34*t35*3.104080344633556E-3;
    double t410 = t45*t46*t47*-3.104080344633556E-3;
    double t411 = t45*t46*t47*3.104080344633556E-3;
    double t412 = t30*t31*t32*-3.105990081579729E-3;
    double t413 = t30*t31*t32*3.105990081579729E-3;
    double t414 = t42*t43*t44*3.105990081579729E-3;
    double t415 = t101+t145;
    double t416 = t104+t146;
    double t417 = t171+t187;
    double t418 = t173+t189;
    double t419 = t85+t246;
    double t420 = t87+t247;
    double t421 = t133+t229;
    double t422 = t135+t231;
    double t423 = t184+t238;
    double t424 = t186+t240;
    double t425 = t66+t112+t159+t214;
    double t426 = t90+t122+t151+t175;
    double t427 = t67+t117+t161+t217;
    double t428 = t93+t123+t154+t179;
    double t429 = t68+t110+t166+t213;
    double t430 = t97+t124+t152+t174;
    double t431 = t69+t118+t169+t216;
    double t432 = t100+t125+t155+t181;
    double t433 = t61+t130+t222+t223;
    double t434 = t81+t82+t196+t243;
    double t435 = t63+t131+t224+t225;
    double t436 = t83+t84+t197+t245;
    double t444 = t263+t363+t374+t382;
    double t445 = t304+t310+t324+t409;
    double t446 = t267+t365+t379+t385;
    double t447 = t308+t314+t328+t410;
    double t448 = t289+t366+t391+t400;
    double t449 = t338+t346+t352+t403;
    double t450 = t293+t370+t395+t401;
    double t451 = t341+t349+t355+t405;
    double t452 = t266+t272+t281+t384;
    double t453 = t251+t312+t327+t333;
    double t454 = t270+t276+t283+t388;
    double t455 = t252+t316+t330+t336;
    double t456 = t286+t292+t298+t368;
    double t457 = t257+t318+t340+t354;
    double t458 = t287+t295+t301+t372;
    double t459 = t259+t322+t344+t357;
    double t460 = t72+t193+t296+t367+t390+t397;
    double t461 = t107+t142+t317+t345+t353+t412;
    double t462 = t74+t195+t300+t371+t394+t398;
    double t463 = t109+t143+t320+t348+t356+t414;
    double t464 = t163+t234+t271+t359+t373+t383;
    double t465 = t203+t210+t303+t311+t331+t406;
    double t466 = t164+t237+t275+t361+t377+t386;
    double t467 = t205+t211+t307+t315+t335+t407;
    double t437 = t433*x13;
    double t438 = t434*x13;
    double t439 = t435*x28;
    double t440 = t436*x28;
    double t441 = -t437;
    double t442 = -t438;
    double t443 = -t439;

    Jx_partial_dq.setZero();
    Jx_partial_dq(0, 28) = -x29*(t18*1.034589188110661E-3+t42*5.699060997402858E-2-t108-t221+t262+t323+t343+t350)-t459*x30+t463*x31;
    Jx_partial_dq(0, 29) = -x30*(t262+t323+t343+t350)+t451*x31-t459*x29;
    Jx_partial_dq(0, 30) = x31*(t108+t221-t262-t323-t343+t351)+t451*x30+t463*x29;
    Jx_partial_dq(0, 34) = t428*x36-x35*(t51+t60-t78-t92+t182+t200);
    Jx_partial_dq(0, 35) = t428*x35+x36*(t78+t92-t182-t200);
    Jx_partial_dq(1, 28) = t458*x30+t462*x31+x29*(t18*-5.699060997402858E-2+t42*1.034589188110661E-3+t149+t194+t279+t294+t302+t396);
    Jx_partial_dq(1, 29) = x30*(t279+t294+t302+t396)+t450*x31+t458*x29;
    Jx_partial_dq(1, 30) = t450*x30+t462*x29+x31*(t149+t194+t279+t294+t302+t396);
    Jx_partial_dq(1, 34) = -t427*x36+x35*(t53+t58+t121+t128+t137+t160);
    Jx_partial_dq(1, 35) = x36*(t121+t128+t137+t160)-t427*x35;
    Jx_partial_dq(2, 29) = -x30*(t80+t190)+t418*x31;
    Jx_partial_dq(2, 30) = -x31*(t44*3.39756885202024E-1+t80+t190)+t418*x30;
    Jx_partial_dq(2, 35) = -x36*(t25*1.79E-2-t56);
    Jx_partial_dq(3, 31) = t455*x33+t467*x34+x32*(t21*-1.034589188110661E-3+t45*5.699060997402856E-2+t204+t250+t255+t309+t329+t337);
    Jx_partial_dq(3, 32) = x33*(t255+t309+t329+t337)+t447*x34+t455*x32;
    Jx_partial_dq(3, 33) = t447*x33+t467*x32+x34*(t204+t250+t255+t309+t329+t337);
    Jx_partial_dq(3, 34) = -t432*x36-x35*(t51+t60-t77+t99+t183-t207);
    Jx_partial_dq(3, 35) = -t432*x35+x36*(t77-t99-t183+t207);
    Jx_partial_dq(4, 31) = -t454*x33-t466*x34-x32*(t21*5.699060997402856E-2+t45*1.034589188110661E-3+t227+t236+t269+t277+t285+t380);
    Jx_partial_dq(4, 32) = -x33*(t269+t277+t285+t380)+t446*x34-t454*x32;
    Jx_partial_dq(4, 33) = t446*x33-t466*x32-x34*(t227+t236+t269+t277+t285+t380);
    Jx_partial_dq(4, 34) = t431*x36+x35*(t53+t58+t120+t129+t140+t167);
    Jx_partial_dq(4, 35) = x36*(t120+t129+t140+t167)+t431*x35;
    Jx_partial_dq(5, 32) = x33*(t95+t241)-t424*x34;
    Jx_partial_dq(5, 33) = x34*(t47*2.875818595898751E-1+t95+t241)-t424*x33;
    Jx_partial_dq(5, 35) = x36*(t25*1.81E-2+t56);
    Jx_partial_dq(6, 13) = t457*x15+t461*x16+x14*(t6*-1.034589188110661E-3+t30*5.699060997402858E-2+t106+t219+t260+t319+t339+t347);
    Jx_partial_dq(6, 14) = x15*(t260+t319+t339+t347)+t449*x16+t457*x14;
    Jx_partial_dq(6, 15) = t449*x15+t461*x14+x16*(t106+t219+t260+t319+t339+t347);
    Jx_partial_dq(6, 19) = -x20*(t50-t59-t76+t89-t177+t199)+t426*x21;
    Jx_partial_dq(6, 20) = t426*x20+x21*(t76-t89+t177-t199);
    Jx_partial_dq(7, 13) = x14*(t6*5.699060997402858E-2+t30*1.034589188110661E-3-t148-t192+t278+t291+t299+t392)+t456*x15-t460*x16;
    Jx_partial_dq(7, 14) = x15*(t278+t291+t299+t392)-t448*x16+t456*x14;
    Jx_partial_dq(7, 15) = -x16*(t148+t192-t278-t291-t299+t393)-t448*x15-t460*x14;
    Jx_partial_dq(7, 19) = t425*x21-x20*(t52+t57+t113+t126+t136+t157);
    Jx_partial_dq(7, 20) = -x21*(t113+t126+t136+t157)+t425*x20;
    Jx_partial_dq(8, 14) = -x15*(t79+t188)-t417*x16;
    Jx_partial_dq(8, 15) = -x16*(t32*3.39756885202024E-1+t79+t188)-t417*x15;
    Jx_partial_dq(8, 20) = -x21*(t13*1.79E-2+t55);
    Jx_partial_dq(9, 16) = -t453*x18+t465*x19-x17*(t9*1.034589188110661E-3+t33*5.699060997402856E-2+t201+t249+t254+t305+t326+t334);
    Jx_partial_dq(9, 17) = -x18*(t254+t305+t326+t334)+t445*x19-t453*x17;
    Jx_partial_dq(9, 18) = t445*x18+t465*x17-x19*(t201+t249+t254+t305+t326+t334);
    Jx_partial_dq(9, 19) = -t430*x21+x20*(-t50+t59+t75+t96+t178+t206);
    Jx_partial_dq(9, 20) = x21*(t75+t96+t178+t206)-t430*x20;
    Jx_partial_dq(10, 16) = -t452*x18-t464*x19+x17*(t9*5.699060997402856E-2-t33*1.034589188110661E-3+t226+t233+t264+t273+t284+t375);
    Jx_partial_dq(10, 17) = -t444*x19-t452*x17-x18*(t265+t274-t284+t376);
    Jx_partial_dq(10, 18) = -t444*x18-t464*x17+x19*(t226+t233+t264+t273+t284+t375);
    Jx_partial_dq(10, 19) = -t429*x21-x20*(t52+t57+t114+t127+t138+t165);
    Jx_partial_dq(10, 20) = -x21*(t114+t127+t138+t165)-t429*x20;
    Jx_partial_dq(11, 17) = x18*(t94+t239)+t423*x19;
    Jx_partial_dq(11, 18) = x19*(t35*2.875818595898751E-1+t94+t239)+t423*x18;
    Jx_partial_dq(11, 20) = x21*(t13*1.81E-2-t55);
    Jx_partial_dq(12, 24) = t422*x26-x25*(t65-t232);
    Jx_partial_dq(12, 25) = t422*x25+x26*(t39*4.365860704565494E-4-t65+t232);
    Jx_partial_dq(12, 26) = t440+x27*(t16*6.740600000000002E-2-t40*4.95436E-1+t436);
    Jx_partial_dq(12, 27) = t440+t436*x27;
    Jx_partial_dq(13, 24) = -x25*(t88+t209)-t420*x26;
    Jx_partial_dq(13, 25) = -x26*(t39*-3.566153386244494E-2+t88+t209)-t420*x25;
    Jx_partial_dq(13, 26) = t443+x27*(t16*4.95436E-1+t40*6.740600000000002E-2-t63+t132-t224-t225);
    Jx_partial_dq(13, 27) = t443-t435*x27;
    Jx_partial_dq(14, 24) = x25*(t71+t105)+t416*x26;
    Jx_partial_dq(14, 25) = t416*x25+x26*(t39*4.987264424463382E-1+t71+t105);
    Jx_partial_dq(15, 9) = -x10*(t64+t230)-t421*x11;
    Jx_partial_dq(15, 10) = -t421*x10-x11*(t27*-4.365860704565494E-4+t64+t230);
    Jx_partial_dq(15, 11) = t442+x12*(t4*6.740600000000002E-2+t28*4.95436E-1-t81-t82-t196+t244);
    Jx_partial_dq(15, 12) = t442-t434*x12;
    Jx_partial_dq(16, 9) = -t419*x11-x10*(t86-t208);
    Jx_partial_dq(16, 10) = -t419*x10-x11*(t27*3.566153386244494E-2+t86-t208);
    Jx_partial_dq(16, 11) = t441-x12*(t4*4.95436E-1-t28*6.740600000000002E-2+t433);
    Jx_partial_dq(16, 12) = t441-t433*x12;
    Jx_partial_dq(17, 9) = x10*(t70+t102)-t415*x11;
    Jx_partial_dq(17, 10) = -t415*x10+x11*(t27*4.987264424463382E-1+t70+t102);

    fkPtr_->compute(0, 
                    contact_joint_id, 
                    q, 
                    nullptr, 
                    &stance_foot_endT, 
                    2);
    
    fkPtr_->getTranslationHessian(H_translation);
    fkPtr_->getRPYHessian(H_rotation);

    Jx_partial_dq.row(18) = H_translation(0) * x;
    Jx_partial_dq.row(19) = H_translation(1) * x;
    Jx_partial_dq.row(20) = H_translation(2) * x;
    Jx_partial_dq.row(21) = H_rotation(0) * x;
    Jx_partial_dq.row(22) = H_rotation(1) * x;
    Jx_partial_dq.row(23) = H_rotation(2) * x;
}

void DigitDynamicsConstraints::get_JTx_partial_dq(const VecX& q, const VecX& x) {
    assert(x.size() == NUM_DEPENDENT_JOINTS);
    assert(JTx_partial_dq.rows() == modelPtr_->nv);
    assert(JTx_partial_dq.cols() == modelPtr_->nv);

    double x1 = x(0);
    double x2 = x(1);
    double x3 = x(2);
    double x4 = x(3);
    double x5 = x(4);
    double x6 = x(5);
    double x7 = x(6);
    double x8 = x(7);
    double x9 = x(8);
    double x10 = x(9);
    double x11 = x(10);
    double x12 = x(11);
    double x13 = x(12);
    double x14 = x(13);
    double x15 = x(14);
    double x16 = x(15);
    double x17 = x(16);
    double x18 = x(17);
    double x19 = x(18);

    double t2 = cos(q(9));
    double t3 = cos(q(10));
    double t4 = cos(q(11));
    double t5 = cos(q(12));
    double t6 = cos(q(13));
    double t7 = cos(q(14));
    double t8 = cos(q(15));
    double t9 = cos(q(16));
    double t10 = cos(q(17));
    double t11 = cos(q(18));
    double t12 = cos(q(19));
    double t13 = cos(q(20));
    double t14 = cos(q(24));
    double t15 = cos(q(25));
    double t16 = cos(q(26));
    double t17 = cos(q(27));
    double t18 = cos(q(28));
    double t19 = cos(q(29));
    double t20 = cos(q(30));
    double t21 = cos(q(31));
    double t22 = cos(q(32));
    double t23 = cos(q(33));
    double t24 = cos(q(34));
    double t25 = cos(q(35));
    double t26 = sin(q(9));
    double t27 = sin(q(10));
    double t28 = sin(q(11));
    double t29 = sin(q(12));
    double t30 = sin(q(13));
    double t31 = sin(q(14));
    double t32 = sin(q(15));
    double t33 = sin(q(16));
    double t34 = sin(q(17));
    double t35 = sin(q(18));
    double t36 = sin(q(19));
    double t37 = sin(q(20));
    double t38 = sin(q(24));
    double t39 = sin(q(25));
    double t40 = sin(q(26));
    double t41 = sin(q(27));
    double t42 = sin(q(28));
    double t43 = sin(q(29));
    double t44 = sin(q(30));
    double t45 = sin(q(31));
    double t46 = sin(q(32));
    double t47 = sin(q(33));
    double t48 = sin(q(34));
    double t49 = sin(q(35));
    double t50 = t12*1.696216709330505E-2;
    double t51 = t24*1.696216709330505E-2;
    double t52 = t36*1.696216709330505E-2;
    double t53 = t48*-1.696216709330505E-2;
    double t54 = t48*1.696216709330505E-2;
    double t55 = t37*9.551E-3;
    double t56 = t49*9.551E-3;
    double t57 = t12*5.143951577823025E-2;
    double t58 = t24*5.143951577823025E-2;
    double t59 = t36*5.143951577823025E-2;
    double t60 = t48*5.143951577823025E-2;
    double t61 = t4*t5*-6.370345907739961E-5;
    double t62 = t4*t5*6.370345907739961E-5;
    double t63 = t16*t17*6.370345907739961E-5;
    double t64 = t2*t3*6.725214316796237E-2;
    double t65 = t14*t15*6.725214316796237E-2;
    double t66 = t12*t13*1.69996184260823E-2;
    double t67 = t24*t25*1.69996184260823E-2;
    double t68 = t12*t13*1.718955829676478E-2;
    double t69 = t24*t25*1.718955829676478E-2;
    double t70 = t2*t3*3.539606431708408E-2;
    double t71 = t14*t15*3.539606431708408E-2;
    double t72 = t6*t8*-2.246221860400801E-3;
    double t73 = t6*t8*2.246221860400801E-3;
    double t74 = t18*t20*2.246221860400801E-3;
    double t75 = t12*t13*9.070578524442012E-3;
    double t76 = t12*t13*9.070578524442013E-3;
    double t77 = t24*t25*9.070578524442012E-3;
    double t78 = t24*t25*9.070578524442013E-3;
    double t79 = t7*t8*2.360206106172353E-3;
    double t80 = t19*t20*2.360206106172353E-3;
    double t81 = t4*t29*6.370345907739961E-5;
    double t82 = t5*t28*6.370345907739961E-5;
    double t83 = t16*t41*6.370345907739961E-5;
    double t84 = t17*t40*6.370345907739961E-5;
    double t85 = t2*t27*6.711175107535902E-2;
    double t86 = t3*t26*6.711175107535902E-2;
    double t87 = t14*t39*6.711175107535902E-2;
    double t88 = t15*t38*6.711175107535902E-2;
    double t89 = t12*t37*1.69996184260823E-2;
    double t90 = t13*t36*-1.69996184260823E-2;
    double t91 = t13*t36*1.69996184260823E-2;
    double t92 = t24*t49*1.69996184260823E-2;
    double t93 = t25*t48*1.69996184260823E-2;
    double t94 = t10*t11*6.574122182670539E-4;
    double t95 = t22*t23*6.574122182670539E-4;
    double t96 = t12*t37*1.718955829676478E-2;
    double t97 = t13*t36*-1.718955829676478E-2;
    double t98 = t13*t36*1.718955829676478E-2;
    double t99 = t24*t49*1.718955829676478E-2;
    double t100 = t25*t48*1.718955829676478E-2;
    double t101 = t2*t27*4.365115769377904E-3;
    double t102 = t3*t26*-4.365115769377904E-3;
    double t103 = t3*t26*4.365115769377904E-3;
    double t104 = t14*t39*4.365115769377904E-3;
    double t105 = t15*t38*4.365115769377904E-3;
    double t106 = t6*t32*2.246221860400801E-3;
    double t107 = t8*t30*2.246221860400801E-3;
    double t108 = t18*t44*2.246221860400801E-3;
    double t109 = t20*t42*2.246221860400801E-3;
    double t110 = t12*t37*-9.070578524442012E-3;
    double t111 = t12*t37*9.070578524442012E-3;
    double t112 = t12*t37*9.070578524442013E-3;
    double t113 = t13*t36*-9.070578524442013E-3;
    double t114 = t13*t36*-9.070578524442012E-3;
    double t115 = t13*t36*9.070578524442012E-3;
    double t116 = t13*t36*9.070578524442013E-3;
    double t117 = t24*t49*-9.070578524442013E-3;
    double t118 = t24*t49*9.070578524442012E-3;
    double t119 = t24*t49*9.070578524442013E-3;
    double t120 = t25*t48*9.070578524442012E-3;
    double t121 = t25*t48*9.070578524442013E-3;
    double t122 = t12*t13*5.605619802270151E-3;
    double t123 = t24*t25*5.605619802270151E-3;
    double t124 = t12*t13*5.668252425759201E-3;
    double t125 = t24*t25*5.668252425759201E-3;
    double t126 = t12*t13*2.991020934719675E-3;
    double t127 = t12*t13*2.991020934719677E-3;
    double t128 = t24*t25*2.991020934719675E-3;
    double t129 = t24*t25*2.991020934719677E-3;
    double t130 = t28*t29*6.370345907739961E-5;
    double t131 = t40*t41*-6.370345907739961E-5;
    double t132 = t40*t41*6.370345907739961E-5;
    double t133 = t26*t27*-6.725214316796237E-2;
    double t134 = t26*t27*6.725214316796237E-2;
    double t135 = t38*t39*6.725214316796237E-2;
    double t136 = t36*t37*1.69996184260823E-2;
    double t137 = t48*t49*1.69996184260823E-2;
    double t138 = t36*t37*-1.718955829676478E-2;
    double t139 = t36*t37*1.718955829676478E-2;
    double t140 = t48*t49*-1.718955829676478E-2;
    double t141 = t48*t49*1.718955829676478E-2;
    double t142 = t6*t8*1.26575449899492E-2;
    double t143 = t18*t20*-1.26575449899492E-2;
    double t144 = t18*t20*1.26575449899492E-2;
    double t145 = t26*t27*3.539606431708408E-2;
    double t146 = t38*t39*-3.539606431708408E-2;
    double t147 = t38*t39*3.539606431708408E-2;
    double t148 = t30*t32*2.246221860400801E-3;
    double t149 = t42*t44*-2.246221860400801E-3;
    double t150 = t42*t44*2.246221860400801E-3;
    double t151 = t36*t37*-9.070578524442013E-3;
    double t152 = t36*t37*9.070578524442012E-3;
    double t153 = t36*t37*9.070578524442013E-3;
    double t154 = t48*t49*-9.070578524442013E-3;
    double t155 = t48*t49*9.070578524442012E-3;
    double t156 = t48*t49*9.070578524442013E-3;
    double t157 = t12*t37*-5.605619802270151E-3;
    double t158 = t12*t37*5.605619802270151E-3;
    double t159 = t13*t36*5.605619802270151E-3;
    double t160 = t24*t49*5.605619802270151E-3;
    double t161 = t25*t48*-5.605619802270151E-3;
    double t162 = t25*t48*5.605619802270151E-3;
    double t163 = t9*t11*8.234792310892687E-4;
    double t164 = t21*t23*8.234792310892687E-4;
    double t165 = t12*t37*5.668252425759201E-3;
    double t166 = t13*t36*5.668252425759201E-3;
    double t167 = t24*t49*-5.668252425759201E-3;
    double t168 = t24*t49*5.668252425759201E-3;
    double t169 = t25*t48*-5.668252425759201E-3;
    double t170 = t25*t48*5.668252425759201E-3;
    double t171 = t31*t32*-2.360206106172353E-3;
    double t172 = t31*t32*2.360206106172353E-3;
    double t173 = t43*t44*2.360206106172353E-3;
    double t174 = t12*t37*-2.991020934719677E-3;
    double t175 = t12*t37*2.991020934719675E-3;
    double t176 = t12*t37*2.991020934719677E-3;
    double t177 = t13*t36*2.991020934719675E-3;
    double t178 = t13*t36*2.991020934719677E-3;
    double t179 = t24*t49*-2.991020934719675E-3;
    double t180 = t24*t49*2.991020934719675E-3;
    double t181 = t24*t49*2.991020934719677E-3;
    double t182 = t25*t48*2.991020934719675E-3;
    double t183 = t25*t48*2.991020934719677E-3;
    double t184 = t34*t35*-6.574122182670539E-4;
    double t185 = t34*t35*6.574122182670539E-4;
    double t186 = t46*t47*6.574122182670539E-4;
    double t187 = t7*t32*1.263678697118524E-2;
    double t188 = t8*t31*1.263678697118524E-2;
    double t189 = t19*t44*1.263678697118524E-2;
    double t190 = t20*t43*-1.263678697118524E-2;
    double t191 = t20*t43*1.263678697118524E-2;
    double t192 = t6*t32*1.26575449899492E-2;
    double t193 = t8*t30*1.26575449899492E-2;
    double t194 = t18*t44*1.26575449899492E-2;
    double t195 = t20*t42*1.26575449899492E-2;
    double t196 = t4*t5*1.20005624647208E-1;
    double t197 = t16*t17*-1.20005624647208E-1;
    double t198 = t16*t17*1.20005624647208E-1;
    double t199 = t36*t37*5.605619802270151E-3;
    double t200 = t48*t49*5.605619802270151E-3;
    double t201 = t9*t35*-8.234792310892687E-4;
    double t202 = t9*t35*8.234792310892687E-4;
    double t203 = t11*t33*8.234792310892687E-4;
    double t204 = t21*t47*8.234792310892687E-4;
    double t205 = t23*t45*8.234792310892687E-4;
    double t206 = t36*t37*5.668252425759201E-3;
    double t207 = t48*t49*5.668252425759201E-3;
    double t208 = t2*t3*4.941905177865888E-1;
    double t209 = t14*t15*4.941905177865888E-1;
    double t210 = t9*t11*1.549180159908108E-2;
    double t211 = t21*t23*-1.549180159908108E-2;
    double t212 = t21*t23*1.549180159908108E-2;
    double t213 = t36*t37*-2.991020934719677E-3;
    double t214 = t36*t37*2.991020934719675E-3;
    double t215 = t36*t37*2.991020934719677E-3;
    double t216 = t48*t49*-2.991020934719677E-3;
    double t217 = t48*t49*2.991020934719675E-3;
    double t218 = t48*t49*2.991020934719677E-3;
    double t219 = t30*t32*-1.26575449899492E-2;
    double t220 = t30*t32*1.26575449899492E-2;
    double t221 = t42*t44*1.26575449899492E-2;
    double t222 = t4*t29*1.20005624647208E-1;
    double t223 = t5*t28*1.20005624647208E-1;
    double t224 = t16*t41*1.20005624647208E-1;
    double t225 = t17*t40*1.20005624647208E-1;
    double t226 = t33*t35*8.234792310892687E-4;
    double t227 = t45*t47*-8.234792310892687E-4;
    double t228 = t45*t47*8.234792310892687E-4;
    double t229 = t2*t27*4.954563135453205E-1;
    double t230 = t3*t26*4.954563135453205E-1;
    double t231 = t14*t39*4.954563135453205E-1;
    double t232 = t15*t38*4.954563135453205E-1;
    double t233 = t9*t35*1.549180159908108E-2;
    double t234 = t11*t33*-1.549180159908108E-2;
    double t235 = t11*t33*1.549180159908108E-2;
    double t236 = t21*t47*1.549180159908108E-2;
    double t237 = t23*t45*1.549180159908108E-2;
    double t238 = t10*t35*1.549973690114125E-2;
    double t239 = t11*t34*1.549973690114125E-2;
    double t240 = t22*t47*1.549973690114125E-2;
    double t241 = t23*t46*-1.549973690114125E-2;
    double t242 = t23*t46*1.549973690114125E-2;
    double t243 = t28*t29*-1.20005624647208E-1;
    double t244 = t28*t29*1.20005624647208E-1;
    double t245 = t40*t41*1.20005624647208E-1;
    double t246 = t26*t27*4.941905177865888E-1;
    double t247 = t38*t39*-4.941905177865888E-1;
    double t248 = t38*t39*4.941905177865888E-1;
    double t249 = t33*t35*1.549180159908108E-2;
    double t250 = t45*t47*1.549180159908108E-2;
    double t251 = t9*t10*t11*2.87566285868891E-1;
    double t252 = t21*t22*t23*-2.87566285868891E-1;
    double t253 = t21*t22*t23*2.87566285868891E-1;
    double t254 = t9*t10*t11*2.879825211612492E-1;
    double t255 = t21*t22*t23*-2.879825211612492E-1;
    double t256 = t21*t22*t23*2.879825211612492E-1;
    double t257 = t6*t7*t8*-3.397508858570615E-1;
    double t258 = t6*t7*t8*3.397508858570615E-1;
    double t259 = t18*t19*t20*3.397508858570615E-1;
    double t260 = t6*t7*t8*-3.399783924207052E-1;
    double t261 = t6*t7*t8*3.399783924207052E-1;
    double t262 = t18*t19*t20*3.399783924207052E-1;
    double t263 = t9*t10*t35*2.87566285868891E-1;
    double t264 = t9*t11*t34*-2.87566285868891E-1;
    double t265 = t9*t11*t34*2.87566285868891E-1;
    double t266 = t10*t11*t33*2.87566285868891E-1;
    double t267 = t21*t22*t47*-2.87566285868891E-1;
    double t268 = t21*t22*t47*2.87566285868891E-1;
    double t269 = t21*t23*t46*2.87566285868891E-1;
    double t270 = t22*t23*t45*2.87566285868891E-1;
    double t271 = t9*t10*t35*2.879825211612492E-1;
    double t272 = t9*t11*t34*2.879825211612492E-1;
    double t273 = t10*t11*t33*-2.879825211612492E-1;
    double t274 = t10*t11*t33*2.879825211612492E-1;
    double t275 = t21*t22*t47*2.879825211612492E-1;
    double t276 = t21*t23*t46*2.879825211612492E-1;
    double t277 = t22*t23*t45*2.879825211612492E-1;
    double t278 = t6*t7*t8*3.020283789547073E-3;
    double t279 = t18*t19*t20*-3.020283789547073E-3;
    double t280 = t18*t19*t20*3.020283789547073E-3;
    double t281 = t9*t10*t11*-3.064210757541298E-3;
    double t282 = t9*t10*t11*3.064210757541298E-3;
    double t283 = t21*t22*t23*3.064210757541298E-3;
    double t284 = t9*t10*t11*3.104080344633556E-3;
    double t285 = t21*t22*t23*3.104080344633556E-3;
    double t286 = t6*t7*t8*3.105990081579729E-3;
    double t287 = t18*t19*t20*-3.105990081579729E-3;
    double t288 = t18*t19*t20*3.105990081579729E-3;
    double t289 = t6*t7*t32*-3.397508858570615E-1;
    double t290 = t6*t7*t32*3.397508858570615E-1;
    double t291 = t6*t8*t31*3.397508858570615E-1;
    double t292 = t7*t8*t30*3.397508858570615E-1;
    double t293 = t18*t19*t44*3.397508858570615E-1;
    double t294 = t18*t20*t43*3.397508858570615E-1;
    double t295 = t19*t20*t42*3.397508858570615E-1;
    double t296 = t6*t7*t32*-3.399783924207052E-1;
    double t297 = t6*t7*t32*3.399783924207052E-1;
    double t298 = t6*t8*t31*3.399783924207052E-1;
    double t299 = t7*t8*t30*3.399783924207052E-1;
    double t300 = t18*t19*t44*3.399783924207052E-1;
    double t301 = t18*t20*t43*3.399783924207052E-1;
    double t302 = t19*t20*t42*3.399783924207052E-1;
    double t303 = t9*t34*t35*2.87566285868891E-1;
    double t304 = t10*t33*t35*2.87566285868891E-1;
    double t305 = t11*t33*t34*-2.87566285868891E-1;
    double t306 = t11*t33*t34*2.87566285868891E-1;
    double t307 = t21*t46*t47*2.87566285868891E-1;
    double t308 = t22*t45*t47*2.87566285868891E-1;
    double t309 = t23*t45*t46*2.87566285868891E-1;
    double t310 = t9*t34*t35*2.879825211612492E-1;
    double t311 = t10*t33*t35*2.879825211612492E-1;
    double t312 = t11*t33*t34*-2.879825211612492E-1;
    double t313 = t11*t33*t34*2.879825211612492E-1;
    double t314 = t21*t46*t47*2.879825211612492E-1;
    double t315 = t22*t45*t47*2.879825211612492E-1;
    double t316 = t23*t45*t46*2.879825211612492E-1;
    double t317 = t6*t7*t32*3.020283789547073E-3;
    double t318 = t6*t8*t31*3.020283789547073E-3;
    double t319 = t7*t8*t30*3.020283789547073E-3;
    double t320 = t18*t19*t44*-3.020283789547073E-3;
    double t321 = t18*t19*t44*3.020283789547073E-3;
    double t322 = t18*t20*t43*3.020283789547073E-3;
    double t323 = t19*t20*t42*3.020283789547073E-3;
    double t324 = t9*t10*t35*-3.064210757541298E-3;
    double t325 = t9*t10*t35*3.064210757541298E-3;
    double t326 = t9*t11*t34*3.064210757541298E-3;
    double t327 = t10*t11*t33*3.064210757541298E-3;
    double t328 = t21*t22*t47*3.064210757541298E-3;
    double t329 = t21*t23*t46*3.064210757541298E-3;
    double t330 = t22*t23*t45*3.064210757541298E-3;
    double t331 = t9*t10*t35*-3.104080344633556E-3;
    double t332 = t9*t10*t35*3.104080344633556E-3;
    double t333 = t9*t11*t34*3.104080344633556E-3;
    double t334 = t10*t11*t33*3.104080344633556E-3;
    double t335 = t21*t22*t47*3.104080344633556E-3;
    double t336 = t21*t23*t46*3.104080344633556E-3;
    double t337 = t22*t23*t45*3.104080344633556E-3;
    double t338 = t6*t7*t32*3.105990081579729E-3;
    double t339 = t6*t8*t31*3.105990081579729E-3;
    double t340 = t7*t8*t30*3.105990081579729E-3;
    double t341 = t18*t19*t44*-3.105990081579729E-3;
    double t342 = t18*t19*t44*3.105990081579729E-3;
    double t343 = t18*t20*t43*3.105990081579729E-3;
    double t344 = t19*t20*t42*3.105990081579729E-3;
    double t345 = t6*t31*t32*3.397508858570615E-1;
    double t346 = t7*t30*t32*3.397508858570615E-1;
    double t347 = t8*t30*t31*3.397508858570615E-1;
    double t348 = t18*t43*t44*3.397508858570615E-1;
    double t349 = t19*t42*t44*3.397508858570615E-1;
    double t350 = t20*t42*t43*-3.397508858570615E-1;
    double t351 = t20*t42*t43*3.397508858570615E-1;
    double t352 = t6*t31*t32*3.399783924207052E-1;
    double t353 = t7*t30*t32*3.399783924207052E-1;
    double t354 = t8*t30*t31*3.399783924207052E-1;
    double t355 = t18*t43*t44*3.399783924207052E-1;
    double t356 = t19*t42*t44*3.399783924207052E-1;
    double t357 = t20*t42*t43*-3.399783924207052E-1;
    double t358 = t20*t42*t43*3.399783924207052E-1;
    double t359 = t33*t34*t35*-2.87566285868891E-1;
    double t360 = t33*t34*t35*2.87566285868891E-1;
    double t361 = t45*t46*t47*-2.87566285868891E-1;
    double t362 = t45*t46*t47*2.87566285868891E-1;
    double t363 = t33*t34*t35*-2.879825211612492E-1;
    double t364 = t33*t34*t35*2.879825211612492E-1;
    double t365 = t45*t46*t47*2.879825211612492E-1;
    double t366 = t6*t31*t32*3.020283789547073E-3;
    double t367 = t7*t30*t32*3.020283789547073E-3;
    double t368 = t8*t30*t31*-3.020283789547073E-3;
    double t369 = t8*t30*t31*3.020283789547073E-3;
    double t370 = t18*t43*t44*3.020283789547073E-3;
    double t371 = t19*t42*t44*3.020283789547073E-3;
    double t372 = t20*t42*t43*3.020283789547073E-3;
    double t373 = t9*t34*t35*3.064210757541298E-3;
    double t374 = t10*t33*t35*3.064210757541298E-3;
    double t375 = t11*t33*t34*-3.064210757541298E-3;
    double t376 = t11*t33*t34*3.064210757541298E-3;
    double t377 = t21*t46*t47*-3.064210757541298E-3;
    double t378 = t21*t46*t47*3.064210757541298E-3;
    double t379 = t22*t45*t47*3.064210757541298E-3;
    double t380 = t23*t45*t46*-3.064210757541298E-3;
    double t381 = t23*t45*t46*3.064210757541298E-3;
    double t382 = t9*t34*t35*3.104080344633556E-3;
    double t383 = t10*t33*t35*3.104080344633556E-3;
    double t384 = t11*t33*t34*3.104080344633556E-3;
    double t385 = t21*t46*t47*3.104080344633556E-3;
    double t386 = t22*t45*t47*-3.104080344633556E-3;
    double t387 = t22*t45*t47*3.104080344633556E-3;
    double t388 = t23*t45*t46*-3.104080344633556E-3;
    double t389 = t23*t45*t46*3.104080344633556E-3;
    double t390 = t6*t31*t32*3.105990081579729E-3;
    double t391 = t7*t30*t32*3.105990081579729E-3;
    double t392 = t8*t30*t31*-3.105990081579729E-3;
    double t393 = t8*t30*t31*3.105990081579729E-3;
    double t394 = t18*t43*t44*3.105990081579729E-3;
    double t395 = t19*t42*t44*3.105990081579729E-3;
    double t396 = t20*t42*t43*3.105990081579729E-3;
    double t397 = t30*t31*t32*3.397508858570615E-1;
    double t398 = t42*t43*t44*-3.397508858570615E-1;
    double t399 = t42*t43*t44*3.397508858570615E-1;
    double t400 = t30*t31*t32*3.399783924207052E-1;
    double t401 = t42*t43*t44*-3.399783924207052E-1;
    double t402 = t42*t43*t44*3.399783924207052E-1;
    double t403 = t30*t31*t32*-3.020283789547073E-3;
    double t404 = t30*t31*t32*3.020283789547073E-3;
    double t405 = t42*t43*t44*3.020283789547073E-3;
    double t406 = t33*t34*t35*3.064210757541298E-3;
    double t407 = t45*t46*t47*-3.064210757541298E-3;
    double t408 = t45*t46*t47*3.064210757541298E-3;
    double t409 = t33*t34*t35*3.104080344633556E-3;
    double t410 = t45*t46*t47*-3.104080344633556E-3;
    double t411 = t45*t46*t47*3.104080344633556E-3;
    double t412 = t30*t31*t32*-3.105990081579729E-3;
    double t413 = t30*t31*t32*3.105990081579729E-3;
    double t414 = t42*t43*t44*3.105990081579729E-3;
    double t415 = t101+t145;
    double t416 = t104+t146;
    double t417 = t171+t187;
    double t418 = t173+t189;
    double t419 = t85+t246;
    double t420 = t87+t247;
    double t421 = t133+t229;
    double t422 = t135+t231;
    double t423 = t184+t238;
    double t424 = t186+t240;
    double t441 = t66+t112+t159+t214;
    double t442 = t90+t122+t151+t175;
    double t443 = t67+t117+t161+t217;
    double t444 = t93+t123+t154+t179;
    double t445 = t68+t110+t166+t213;
    double t446 = t97+t124+t152+t174;
    double t447 = t69+t118+t169+t216;
    double t448 = t100+t125+t155+t181;
    double t449 = t61+t130+t222+t223;
    double t450 = t81+t82+t196+t243;
    double t451 = t63+t131+t224+t225;
    double t452 = t83+t84+t197+t245;
    double t472 = t263+t363+t374+t382;
    double t473 = t304+t310+t324+t409;
    double t474 = t267+t365+t379+t385;
    double t475 = t308+t314+t328+t410;
    double t476 = t289+t366+t391+t400;
    double t477 = t338+t346+t352+t403;
    double t478 = t293+t370+t395+t401;
    double t479 = t341+t349+t355+t405;
    double t480 = t266+t272+t281+t384;
    double t481 = t251+t312+t327+t333;
    double t482 = t270+t276+t283+t388;
    double t483 = t252+t316+t330+t336;
    double t484 = t286+t292+t298+t368;
    double t485 = t257+t318+t340+t354;
    double t486 = t287+t295+t301+t372;
    double t487 = t259+t322+t344+t357;
    double t510 = t72+t193+t296+t367+t390+t397;
    double t511 = t107+t142+t317+t345+t353+t412;
    double t512 = t74+t195+t300+t371+t394+t398;
    double t513 = t109+t143+t320+t348+t356+t414;
    double t514 = t163+t234+t271+t359+t373+t383;
    double t515 = t203+t210+t303+t311+t331+t406;
    double t516 = t164+t237+t275+t361+t377+t386;
    double t517 = t205+t211+t307+t315+t335+t407;
    double t425 = t415*x18;
    double t426 = t416*x15;
    double t427 = t417*x9;
    double t428 = t418*x3;
    double t429 = t419*x17;
    double t430 = t420*x14;
    double t431 = t421*x16;
    double t432 = t422*x13;
    double t433 = t423*x12;
    double t434 = t424*x6;
    double t453 = t442*x7;
    double t454 = t441*x8;
    double t455 = t444*x1;
    double t456 = t443*x2;
    double t457 = t446*x10;
    double t458 = t445*x11;
    double t459 = t448*x4;
    double t460 = t447*x5;
    double t461 = t450*x16;
    double t462 = t449*x17;
    double t463 = t452*x13;
    double t464 = t451*x14;
    double t488 = t481*x10;
    double t489 = t480*x11;
    double t490 = t483*x4;
    double t491 = t482*x5;
    double t492 = t485*x7;
    double t493 = t484*x8;
    double t494 = t487*x1;
    double t495 = t486*x2;
    double t496 = t473*x10;
    double t497 = t472*x11;
    double t498 = t475*x4;
    double t499 = t474*x5;
    double t500 = t477*x7;
    double t501 = t476*x8;
    double t502 = t479*x1;
    double t503 = t478*x2;
    double t518 = t515*x10;
    double t519 = t514*x11;
    double t520 = t517*x4;
    double t521 = t516*x5;
    double t522 = t511*x7;
    double t523 = t510*x8;
    double t524 = t513*x1;
    double t525 = t512*x2;
    double t435 = -t425;
    double t436 = -t427;
    double t437 = -t429;
    double t438 = -t430;
    double t439 = -t431;
    double t440 = -t434;
    double t465 = -t456;
    double t466 = -t457;
    double t467 = -t458;
    double t468 = -t459;
    double t469 = -t461;
    double t470 = -t462;
    double t471 = -t464;
    double t504 = -t488;
    double t505 = -t489;
    double t506 = -t491;
    double t507 = -t494;
    double t508 = -t497;
    double t509 = -t501;
    double t527 = -t519;
    double t528 = -t521;
    double t529 = -t523;
    double t533 = t492+t493;
    double t537 = t428+t502+t503;
    double t541 = t524+t525;
    double t526 = t426+t432+t438;
    double t530 = t435+t437+t439;
    double t531 = t463+t471;
    double t532 = t469+t470;
    double t534 = t490+t506;
    double t535 = t495+t507;
    double t536 = t504+t505;
    double t538 = t433+t496+t508;
    double t539 = t440+t498+t499;
    double t540 = t436+t500+t509;
    double t542 = t522+t529;
    double t543 = t518+t527;
    double t544 = t520+t528;
    double t545 = t453+t454+t466+t467;
    double t546 = t455+t460+t465+t468;

    JTx_partial_dq.setZero();
    JTx_partial_dq(9, 9) = x18*(t70+t102)-x16*(t64+t230)-x17*(t86-t208);
    JTx_partial_dq(9, 10) = t530;
    JTx_partial_dq(10, 9) = t530;
    JTx_partial_dq(10, 10) = -x16*(t27*-4.365860704565494E-4+t64+t230)+x18*(t27*4.987264424463382E-1+t70+t102)-x17*(t27*3.566153386244494E-2+t86-t208);
    JTx_partial_dq(11, 11) = -x17*(t4*4.95436E-1-t28*6.740600000000002E-2+t449)+x16*(t4*6.740600000000002E-2+t28*4.95436E-1-t81-t82-t196+t244);
    JTx_partial_dq(11, 12) = t532;
    JTx_partial_dq(12, 11) = t532;
    JTx_partial_dq(12, 12) = t532;
    JTx_partial_dq(13, 13) = x8*(t6*5.699060997402858E-2+t30*1.034589188110661E-3-t148-t192+t278+t291+t299+t392)+x7*(t6*-1.034589188110661E-3+t30*5.699060997402858E-2+t106+t219+t260+t319+t339+t347);
    JTx_partial_dq(13, 14) = t533;
    JTx_partial_dq(13, 15) = t542;
    JTx_partial_dq(14, 13) = t533;
    JTx_partial_dq(14, 14) = x8*(t278+t291+t299+t392)+x7*(t260+t319+t339+t347)-x9*(t79+t188);
    JTx_partial_dq(14, 15) = t540;
    JTx_partial_dq(15, 13) = t542;
    JTx_partial_dq(15, 14) = t540;
    JTx_partial_dq(15, 15) = -x8*(t148+t192-t278-t291-t299+t393)-x9*(t32*3.39756885202024E-1+t79+t188)+x7*(t106+t219+t260+t319+t339+t347);
    JTx_partial_dq(16, 16) = x11*(t9*5.699060997402856E-2-t33*1.034589188110661E-3+t226+t233+t264+t273+t284+t375)-x10*(t9*1.034589188110661E-3+t33*5.699060997402856E-2+t201+t249+t254+t305+t326+t334);
    JTx_partial_dq(16, 17) = t536;
    JTx_partial_dq(16, 18) = t543;
    JTx_partial_dq(17, 16) = t536;
    JTx_partial_dq(17, 17) = -x10*(t254+t305+t326+t334)+x12*(t94+t239)-x11*(t265+t274-t284+t376);
    JTx_partial_dq(17, 18) = t538;
    JTx_partial_dq(18, 16) = t543;
    JTx_partial_dq(18, 17) = t538;
    JTx_partial_dq(18, 18) = x12*(t35*2.875818595898751E-1+t94+t239)+x11*(t226+t233+t264+t273+t284+t375)-x10*(t201+t249+t254+t305+t326+t334);
    JTx_partial_dq(19, 19) = -x7*(t50-t59-t76+t89-t177+t199)+x10*(-t50+t59+t75+t96+t178+t206)-x8*(t52+t57+t113+t126+t136+t157)-x11*(t52+t57+t114+t127+t138+t165);
    JTx_partial_dq(19, 20) = t545;
    JTx_partial_dq(20, 19) = t545;
    JTx_partial_dq(20, 20) = -x8*(t113+t126+t136+t157)-x11*(t114+t127+t138+t165)+x10*(t75+t96+t178+t206)+x12*(t13*1.81E-2-t55)+x7*(t76-t89+t177-t199)-x9*(t13*1.79E-2+t55);
    JTx_partial_dq(24, 24) = x15*(t71+t105)-x14*(t88+t209)-x13*(t65-t232);
    JTx_partial_dq(24, 25) = t526;
    JTx_partial_dq(25, 24) = t526;
    JTx_partial_dq(25, 25) = -x14*(t39*-3.566153386244494E-2+t88+t209)+x15*(t39*4.987264424463382E-1+t71+t105)+x13*(t39*4.365860704565494E-4-t65+t232);
    JTx_partial_dq(26, 26) = x13*(t16*6.740600000000002E-2-t40*4.95436E-1+t452)+x14*(t16*4.95436E-1+t40*6.740600000000002E-2-t63+t132-t224-t225);
    JTx_partial_dq(26, 27) = t531;
    JTx_partial_dq(27, 26) = t531;
    JTx_partial_dq(27, 27) = t531;
    JTx_partial_dq(28, 28) = -x1*(t18*1.034589188110661E-3+t42*5.699060997402858E-2-t108-t221+t262+t323+t343+t350)+x2*(t18*-5.699060997402858E-2+t42*1.034589188110661E-3+t149+t194+t279+t294+t302+t396);
    JTx_partial_dq(28, 29) = t535;
    JTx_partial_dq(28, 30) = t541;
    JTx_partial_dq(29, 28) = t535;
    JTx_partial_dq(29, 29) = x2*(t279+t294+t302+t396)-x1*(t262+t323+t343+t350)-x3*(t80+t190);
    JTx_partial_dq(29, 30) = t537;
    JTx_partial_dq(30, 28) = t541;
    JTx_partial_dq(30, 29) = t537;
    JTx_partial_dq(30, 30) = x1*(t108+t221-t262-t323-t343+t351)-x3*(t44*3.39756885202024E-1+t80+t190)+x2*(t149+t194+t279+t294+t302+t396);
    JTx_partial_dq(31, 31) = -x5*(t21*5.699060997402856E-2+t45*1.034589188110661E-3+t227+t236+t269+t277+t285+t380)+x4*(t21*-1.034589188110661E-3+t45*5.699060997402856E-2+t204+t250+t255+t309+t329+t337);
    JTx_partial_dq(31, 32) = t534;
    JTx_partial_dq(31, 33) = t544;
    JTx_partial_dq(32, 31) = t534;
    JTx_partial_dq(32, 32) = -x5*(t269+t277+t285+t380)+x4*(t255+t309+t329+t337)+x6*(t95+t241);
    JTx_partial_dq(32, 33) = t539;
    JTx_partial_dq(33, 31) = t544;
    JTx_partial_dq(33, 32) = t539;
    JTx_partial_dq(33, 33) = x6*(t47*2.875818595898751E-1+t95+t241)-x5*(t227+t236+t269+t277+t285+t380)+x4*(t204+t250+t255+t309+t329+t337);
    JTx_partial_dq(34, 34) = -x1*(t51+t60-t78-t92+t182+t200)-x4*(t51+t60-t77+t99+t183-t207)+x2*(t53+t58+t121+t128+t137+t160)+x5*(t53+t58+t120+t129+t140+t167);
    JTx_partial_dq(34, 35) = t546;
    JTx_partial_dq(35, 34) = t546;
    JTx_partial_dq(35, 35) = x2*(t121+t128+t137+t160)+x5*(t120+t129+t140+t167)-x3*(t25*1.79E-2-t56)+x1*(t78+t92-t182-t200)+x4*(t77-t99-t183+t207)+x6*(t25*1.81E-2+t56);

    fkPtr_->compute(0, 
                    contact_joint_id, 
                    q, 
                    nullptr, 
                    &stance_foot_endT, 
                    2);

    fkPtr_->getTranslationHessian(H_translation);
    fkPtr_->getRPYHessian(H_rotation);

    JTx_partial_dq += H_translation(0) * x(18);
    JTx_partial_dq += H_translation(1) * x(19);
    JTx_partial_dq += H_translation(2) * x(20);
    JTx_partial_dq += H_rotation(0) * x(21);
    JTx_partial_dq += H_rotation(1) * x(22);
    JTx_partial_dq += H_rotation(2) * x(23);
}

void DigitDynamicsConstraints::get_Jxy_partial_dq(const VecX& q, const VecX& x, const VecX& y) {
    assert(x.size() == modelPtr_->nv);
    assert(y.size() == modelPtr_->nv);
    assert(Jxy_partial_dq.rows() == NUM_DEPENDENT_JOINTS);
    assert(Jxy_partial_dq.cols() == modelPtr_->nv);

    double x1 = x(0);
    double x2 = x(1);
    double x3 = x(2);
    double x4 = x(3);
    double x5 = x(4);
    double x6 = x(5);
    double x7 = x(6);
    double x8 = x(7);
    double x9 = x(8);
    double x10 = x(9);
    double x11 = x(10);
    double x12 = x(11);
    double x13 = x(12);
    double x14 = x(13);
    double x15 = x(14);
    double x16 = x(15);
    double x17 = x(16);
    double x18 = x(17);
    double x19 = x(18);
    double x20 = x(19);
    double x21 = x(20);
    double x22 = x(21);
    double x23 = x(22);
    double x24 = x(23);
    double x25 = x(24);
    double x26 = x(25);
    double x27 = x(26);
    double x28 = x(27);
    double x29 = x(28);
    double x30 = x(29);
    double x31 = x(30);
    double x32 = x(31);
    double x33 = x(32);
    double x34 = x(33);
    double x35 = x(34);
    double x36 = x(35);

    double y1 = y(0);
    double y2 = y(1);
    double y3 = y(2);
    double y4 = y(3);
    double y5 = y(4);
    double y6 = y(5);
    double y7 = y(6);
    double y8 = y(7);
    double y9 = y(8);
    double y10 = y(9);
    double y11 = y(10);
    double y12 = y(11);
    double y13 = y(12);
    double y14 = y(13);
    double y15 = y(14);
    double y16 = y(15);
    double y17 = y(16);
    double y18 = y(17);
    double y19 = y(18);
    double y20 = y(19);
    double y21 = y(20);
    double y22 = y(21);
    double y23 = y(22);
    double y24 = y(23);
    double y25 = y(24);
    double y26 = y(25);
    double y27 = y(26);
    double y28 = y(27);
    double y29 = y(28);
    double y30 = y(29);
    double y31 = y(30);
    double y32 = y(31);
    double y33 = y(32);
    double y34 = y(33);
    double y35 = y(34);
    double y36 = y(35);

    double t2 = cos(q(9));
    double t3 = cos(q(10));
    double t4 = cos(q(11));
    double t5 = cos(q(12));
    double t6 = cos(q(13));
    double t7 = cos(q(14));
    double t8 = cos(q(15));
    double t9 = cos(q(16));
    double t10 = cos(q(17));
    double t11 = cos(q(18));
    double t12 = cos(q(19));
    double t13 = cos(q(20));
    double t14 = cos(q(24));
    double t15 = cos(q(25));
    double t16 = cos(q(26));
    double t17 = cos(q(27));
    double t18 = cos(q(28));
    double t19 = cos(q(29));
    double t20 = cos(q(30));
    double t21 = cos(q(31));
    double t22 = cos(q(32));
    double t23 = cos(q(33));
    double t24 = cos(q(34));
    double t25 = cos(q(35));
    double t26 = sin(q(9));
    double t27 = sin(q(10));
    double t28 = sin(q(11));
    double t29 = sin(q(12));
    double t30 = sin(q(13));
    double t31 = sin(q(14));
    double t32 = sin(q(15));
    double t33 = sin(q(16));
    double t34 = sin(q(17));
    double t35 = sin(q(18));
    double t36 = sin(q(19));
    double t37 = sin(q(20));
    double t38 = sin(q(24));
    double t39 = sin(q(25));
    double t40 = sin(q(26));
    double t41 = sin(q(27));
    double t42 = sin(q(28));
    double t43 = sin(q(29));
    double t44 = sin(q(30));
    double t45 = sin(q(31));
    double t46 = sin(q(32));
    double t47 = sin(q(33));
    double t48 = sin(q(34));
    double t49 = sin(q(35));
    double t50 = t12*1.696216709330505E-2;
    double t51 = t24*1.696216709330505E-2;
    double t52 = t13*9.551E-3;
    double t53 = t25*9.551E-3;
    double t54 = t36*1.696216709330505E-2;
    double t55 = t48*-1.696216709330505E-2;
    double t56 = t48*1.696216709330505E-2;
    double t57 = t12*5.143951577823025E-2;
    double t58 = t24*5.143951577823025E-2;
    double t59 = t36*5.143951577823025E-2;
    double t60 = t48*5.143951577823025E-2;
    double t61 = t4*t5*-6.370345907739961E-5;
    double t62 = t4*t5*6.370345907739961E-5;
    double t63 = t16*t17*6.370345907739961E-5;
    double t64 = t2*t3*6.711175107535902E-2;
    double t65 = t14*t15*6.711175107535902E-2;
    double t66 = t12*t13*1.69996184260823E-2;
    double t67 = t24*t25*1.69996184260823E-2;
    double t68 = t12*t13*1.718955829676478E-2;
    double t69 = t24*t25*1.718955829676478E-2;
    double t70 = t2*t3*4.365115769377904E-3;
    double t71 = t14*t15*4.365115769377904E-3;
    double t72 = t6*t8*-2.246221860400801E-3;
    double t73 = t6*t8*2.246221860400801E-3;
    double t74 = t18*t20*2.246221860400801E-3;
    double t75 = t12*t13*9.070578524442012E-3;
    double t76 = t12*t13*9.070578524442013E-3;
    double t77 = t24*t25*9.070578524442012E-3;
    double t78 = t24*t25*9.070578524442013E-3;
    double t79 = t4*t29*6.370345907739961E-5;
    double t80 = t5*t28*6.370345907739961E-5;
    double t81 = t16*t41*6.370345907739961E-5;
    double t82 = t17*t40*6.370345907739961E-5;
    double t83 = t2*t27*6.725214316796237E-2;
    double t84 = t3*t26*-6.725214316796237E-2;
    double t85 = t3*t26*6.725214316796237E-2;
    double t86 = t14*t39*6.725214316796237E-2;
    double t87 = t15*t38*6.725214316796237E-2;
    double t88 = t12*t37*-1.69996184260823E-2;
    double t89 = t12*t37*1.69996184260823E-2;
    double t90 = t13*t36*-1.69996184260823E-2;
    double t91 = t13*t36*1.69996184260823E-2;
    double t92 = t24*t49*1.69996184260823E-2;
    double t93 = t25*t48*1.69996184260823E-2;
    double t94 = t12*t37*1.718955829676478E-2;
    double t95 = t13*t36*-1.718955829676478E-2;
    double t96 = t13*t36*1.718955829676478E-2;
    double t97 = t24*t49*-1.718955829676478E-2;
    double t98 = t24*t49*1.718955829676478E-2;
    double t99 = t25*t48*1.718955829676478E-2;
    double t100 = t2*t27*3.539606431708408E-2;
    double t101 = t3*t26*3.539606431708408E-2;
    double t102 = t14*t39*3.539606431708408E-2;
    double t103 = t15*t38*-3.539606431708408E-2;
    double t104 = t15*t38*3.539606431708408E-2;
    double t105 = t6*t32*2.246221860400801E-3;
    double t106 = t8*t30*2.246221860400801E-3;
    double t107 = t18*t44*2.246221860400801E-3;
    double t108 = t20*t42*2.246221860400801E-3;
    double t109 = t12*t37*-9.070578524442012E-3;
    double t110 = t12*t37*9.070578524442012E-3;
    double t111 = t12*t37*9.070578524442013E-3;
    double t112 = t13*t36*-9.070578524442013E-3;
    double t113 = t13*t36*-9.070578524442012E-3;
    double t114 = t13*t36*9.070578524442012E-3;
    double t115 = t13*t36*9.070578524442013E-3;
    double t116 = t24*t49*-9.070578524442013E-3;
    double t117 = t24*t49*9.070578524442012E-3;
    double t118 = t24*t49*9.070578524442013E-3;
    double t119 = t25*t48*9.070578524442012E-3;
    double t120 = t25*t48*9.070578524442013E-3;
    double t121 = t12*t13*5.605619802270151E-3;
    double t122 = t24*t25*5.605619802270151E-3;
    double t123 = t12*t13*5.668252425759201E-3;
    double t124 = t24*t25*5.668252425759201E-3;
    double t125 = t7*t32*2.360206106172353E-3;
    double t126 = t8*t31*-2.360206106172353E-3;
    double t127 = t8*t31*2.360206106172353E-3;
    double t128 = t19*t44*2.360206106172353E-3;
    double t129 = t20*t43*2.360206106172353E-3;
    double t130 = t12*t13*2.991020934719675E-3;
    double t131 = t12*t13*2.991020934719677E-3;
    double t132 = t24*t25*2.991020934719675E-3;
    double t133 = t24*t25*2.991020934719677E-3;
    double t134 = t28*t29*6.370345907739961E-5;
    double t135 = t40*t41*-6.370345907739961E-5;
    double t136 = t40*t41*6.370345907739961E-5;
    double t137 = t26*t27*-6.711175107535902E-2;
    double t138 = t26*t27*6.711175107535902E-2;
    double t139 = t38*t39*6.711175107535902E-2;
    double t140 = t36*t37*1.69996184260823E-2;
    double t141 = t48*t49*1.69996184260823E-2;
    double t142 = t10*t35*6.574122182670539E-4;
    double t143 = t11*t34*-6.574122182670539E-4;
    double t144 = t11*t34*6.574122182670539E-4;
    double t145 = t22*t47*6.574122182670539E-4;
    double t146 = t23*t46*6.574122182670539E-4;
    double t147 = t36*t37*-1.718955829676478E-2;
    double t148 = t36*t37*1.718955829676478E-2;
    double t149 = t48*t49*-1.718955829676478E-2;
    double t150 = t48*t49*1.718955829676478E-2;
    double t151 = t26*t27*-4.365115769377904E-3;
    double t152 = t26*t27*4.365115769377904E-3;
    double t153 = t38*t39*4.365115769377904E-3;
    double t154 = t7*t8*1.263678697118524E-2;
    double t155 = t19*t20*1.263678697118524E-2;
    double t156 = t6*t8*1.26575449899492E-2;
    double t157 = t18*t20*-1.26575449899492E-2;
    double t158 = t18*t20*1.26575449899492E-2;
    double t159 = t30*t32*2.246221860400801E-3;
    double t160 = t42*t44*-2.246221860400801E-3;
    double t161 = t42*t44*2.246221860400801E-3;
    double t162 = t36*t37*-9.070578524442013E-3;
    double t163 = t36*t37*9.070578524442012E-3;
    double t164 = t36*t37*9.070578524442013E-3;
    double t165 = t48*t49*-9.070578524442013E-3;
    double t166 = t48*t49*9.070578524442012E-3;
    double t167 = t48*t49*9.070578524442013E-3;
    double t168 = t12*t37*-5.605619802270151E-3;
    double t169 = t12*t37*5.605619802270151E-3;
    double t170 = t13*t36*5.605619802270151E-3;
    double t171 = t24*t49*5.605619802270151E-3;
    double t172 = t25*t48*-5.605619802270151E-3;
    double t173 = t25*t48*5.605619802270151E-3;
    double t174 = t9*t11*8.234792310892687E-4;
    double t175 = t21*t23*8.234792310892687E-4;
    double t176 = t12*t37*5.668252425759201E-3;
    double t177 = t13*t36*5.668252425759201E-3;
    double t178 = t24*t49*-5.668252425759201E-3;
    double t179 = t24*t49*5.668252425759201E-3;
    double t180 = t25*t48*-5.668252425759201E-3;
    double t181 = t25*t48*5.668252425759201E-3;
    double t182 = t12*t37*-2.991020934719677E-3;
    double t183 = t12*t37*2.991020934719675E-3;
    double t184 = t12*t37*2.991020934719677E-3;
    double t185 = t13*t36*2.991020934719675E-3;
    double t186 = t13*t36*2.991020934719677E-3;
    double t187 = t24*t49*-2.991020934719675E-3;
    double t188 = t24*t49*2.991020934719675E-3;
    double t189 = t24*t49*2.991020934719677E-3;
    double t190 = t25*t48*-2.991020934719677E-3;
    double t191 = t25*t48*-2.991020934719675E-3;
    double t192 = t25*t48*2.991020934719675E-3;
    double t193 = t25*t48*2.991020934719677E-3;
    double t194 = t6*t32*1.26575449899492E-2;
    double t195 = t8*t30*1.26575449899492E-2;
    double t196 = t18*t44*1.26575449899492E-2;
    double t197 = t20*t42*1.26575449899492E-2;
    double t198 = t4*t5*1.20005624647208E-1;
    double t199 = t16*t17*-1.20005624647208E-1;
    double t200 = t16*t17*1.20005624647208E-1;
    double t201 = t36*t37*-5.605619802270151E-3;
    double t202 = t36*t37*5.605619802270151E-3;
    double t203 = t48*t49*-5.605619802270151E-3;
    double t204 = t48*t49*5.605619802270151E-3;
    double t205 = t9*t35*-8.234792310892687E-4;
    double t206 = t9*t35*8.234792310892687E-4;
    double t207 = t11*t33*8.234792310892687E-4;
    double t208 = t21*t47*8.234792310892687E-4;
    double t209 = t23*t45*8.234792310892687E-4;
    double t210 = t36*t37*5.668252425759201E-3;
    double t211 = t48*t49*5.668252425759201E-3;
    double t212 = t2*t3*4.954563135453205E-1;
    double t213 = t14*t15*4.954563135453205E-1;
    double t214 = t9*t11*1.549180159908108E-2;
    double t215 = t21*t23*-1.549180159908108E-2;
    double t216 = t21*t23*1.549180159908108E-2;
    double t217 = t10*t11*1.549973690114125E-2;
    double t218 = t22*t23*1.549973690114125E-2;
    double t219 = t36*t37*-2.991020934719677E-3;
    double t220 = t36*t37*2.991020934719675E-3;
    double t221 = t36*t37*2.991020934719677E-3;
    double t222 = t48*t49*-2.991020934719677E-3;
    double t223 = t48*t49*2.991020934719675E-3;
    double t224 = t48*t49*2.991020934719677E-3;
    double t225 = t31*t32*1.263678697118524E-2;
    double t226 = t43*t44*-1.263678697118524E-2;
    double t227 = t43*t44*1.263678697118524E-2;
    double t228 = t30*t32*-1.26575449899492E-2;
    double t229 = t30*t32*1.26575449899492E-2;
    double t230 = t42*t44*1.26575449899492E-2;
    double t231 = t4*t29*1.20005624647208E-1;
    double t232 = t5*t28*1.20005624647208E-1;
    double t233 = t16*t41*1.20005624647208E-1;
    double t234 = t17*t40*1.20005624647208E-1;
    double t235 = t33*t35*8.234792310892687E-4;
    double t236 = t45*t47*-8.234792310892687E-4;
    double t237 = t45*t47*8.234792310892687E-4;
    double t238 = t2*t27*4.941905177865888E-1;
    double t239 = t3*t26*4.941905177865888E-1;
    double t240 = t14*t39*4.941905177865888E-1;
    double t241 = t15*t38*-4.941905177865888E-1;
    double t242 = t15*t38*4.941905177865888E-1;
    double t243 = t9*t35*1.549180159908108E-2;
    double t244 = t11*t33*-1.549180159908108E-2;
    double t245 = t11*t33*1.549180159908108E-2;
    double t246 = t21*t47*1.549180159908108E-2;
    double t247 = t23*t45*1.549180159908108E-2;
    double t248 = t28*t29*-1.20005624647208E-1;
    double t249 = t28*t29*1.20005624647208E-1;
    double t250 = t40*t41*1.20005624647208E-1;
    double t251 = t26*t27*4.954563135453205E-1;
    double t252 = t38*t39*-4.954563135453205E-1;
    double t253 = t38*t39*4.954563135453205E-1;
    double t254 = t33*t35*1.549180159908108E-2;
    double t255 = t45*t47*1.549180159908108E-2;
    double t256 = t34*t35*1.549973690114125E-2;
    double t257 = t46*t47*-1.549973690114125E-2;
    double t258 = t46*t47*1.549973690114125E-2;
    double t259 = t9*t10*t11*2.87566285868891E-1;
    double t260 = t21*t22*t23*-2.87566285868891E-1;
    double t261 = t21*t22*t23*2.87566285868891E-1;
    double t262 = t9*t10*t11*2.879825211612492E-1;
    double t263 = t21*t22*t23*-2.879825211612492E-1;
    double t264 = t21*t22*t23*2.879825211612492E-1;
    double t265 = t6*t7*t8*-3.397508858570615E-1;
    double t266 = t6*t7*t8*3.397508858570615E-1;
    double t267 = t18*t19*t20*3.397508858570615E-1;
    double t268 = t6*t7*t8*-3.399783924207052E-1;
    double t269 = t6*t7*t8*3.399783924207052E-1;
    double t270 = t18*t19*t20*-3.399783924207052E-1;
    double t271 = t18*t19*t20*3.399783924207052E-1;
    double t272 = t9*t10*t35*2.87566285868891E-1;
    double t273 = t9*t11*t34*-2.87566285868891E-1;
    double t274 = t9*t11*t34*2.87566285868891E-1;
    double t275 = t10*t11*t33*2.87566285868891E-1;
    double t276 = t21*t22*t47*-2.87566285868891E-1;
    double t277 = t21*t22*t47*2.87566285868891E-1;
    double t278 = t21*t23*t46*2.87566285868891E-1;
    double t279 = t22*t23*t45*2.87566285868891E-1;
    double t280 = t9*t10*t35*2.879825211612492E-1;
    double t281 = t9*t11*t34*2.879825211612492E-1;
    double t282 = t10*t11*t33*-2.879825211612492E-1;
    double t283 = t10*t11*t33*2.879825211612492E-1;
    double t284 = t21*t22*t47*-2.879825211612492E-1;
    double t285 = t21*t22*t47*2.879825211612492E-1;
    double t286 = t21*t23*t46*2.879825211612492E-1;
    double t287 = t22*t23*t45*2.879825211612492E-1;
    double t288 = t6*t7*t8*-3.020283789547073E-3;
    double t289 = t6*t7*t8*3.020283789547073E-3;
    double t290 = t18*t19*t20*-3.020283789547073E-3;
    double t291 = t18*t19*t20*3.020283789547073E-3;
    double t292 = t9*t10*t11*-3.064210757541298E-3;
    double t293 = t9*t10*t11*3.064210757541298E-3;
    double t294 = t21*t22*t23*3.064210757541298E-3;
    double t295 = t9*t10*t11*-3.104080344633556E-3;
    double t296 = t9*t10*t11*3.104080344633556E-3;
    double t297 = t21*t22*t23*3.104080344633556E-3;
    double t298 = t6*t7*t8*3.105990081579729E-3;
    double t299 = t18*t19*t20*-3.105990081579729E-3;
    double t300 = t18*t19*t20*3.105990081579729E-3;
    double t301 = t6*t7*t32*-3.397508858570615E-1;
    double t302 = t6*t7*t32*3.397508858570615E-1;
    double t303 = t6*t8*t31*-3.397508858570615E-1;
    double t304 = t6*t8*t31*3.397508858570615E-1;
    double t305 = t7*t8*t30*3.397508858570615E-1;
    double t306 = t18*t19*t44*3.397508858570615E-1;
    double t307 = t18*t20*t43*3.397508858570615E-1;
    double t308 = t19*t20*t42*3.397508858570615E-1;
    double t309 = t6*t7*t32*-3.399783924207052E-1;
    double t310 = t6*t7*t32*3.399783924207052E-1;
    double t311 = t6*t8*t31*3.399783924207052E-1;
    double t312 = t7*t8*t30*-3.399783924207052E-1;
    double t313 = t7*t8*t30*3.399783924207052E-1;
    double t314 = t18*t19*t44*3.399783924207052E-1;
    double t315 = t18*t20*t43*3.399783924207052E-1;
    double t316 = t19*t20*t42*3.399783924207052E-1;
    double t317 = t9*t34*t35*2.87566285868891E-1;
    double t318 = t10*t33*t35*2.87566285868891E-1;
    double t319 = t11*t33*t34*-2.87566285868891E-1;
    double t320 = t11*t33*t34*2.87566285868891E-1;
    double t321 = t21*t46*t47*2.87566285868891E-1;
    double t322 = t22*t45*t47*2.87566285868891E-1;
    double t323 = t23*t45*t46*2.87566285868891E-1;
    double t324 = t9*t34*t35*2.879825211612492E-1;
    double t325 = t10*t33*t35*2.879825211612492E-1;
    double t326 = t11*t33*t34*-2.879825211612492E-1;
    double t327 = t11*t33*t34*2.879825211612492E-1;
    double t328 = t21*t46*t47*2.879825211612492E-1;
    double t329 = t22*t45*t47*2.879825211612492E-1;
    double t330 = t23*t45*t46*2.879825211612492E-1;
    double t331 = t6*t7*t32*3.020283789547073E-3;
    double t332 = t6*t8*t31*3.020283789547073E-3;
    double t333 = t7*t8*t30*3.020283789547073E-3;
    double t334 = t18*t19*t44*-3.020283789547073E-3;
    double t335 = t18*t19*t44*3.020283789547073E-3;
    double t336 = t18*t20*t43*3.020283789547073E-3;
    double t337 = t19*t20*t42*-3.020283789547073E-3;
    double t338 = t19*t20*t42*3.020283789547073E-3;
    double t339 = t9*t10*t35*-3.064210757541298E-3;
    double t340 = t9*t10*t35*3.064210757541298E-3;
    double t341 = t9*t11*t34*3.064210757541298E-3;
    double t342 = t10*t11*t33*3.064210757541298E-3;
    double t343 = t21*t22*t47*3.064210757541298E-3;
    double t344 = t21*t23*t46*3.064210757541298E-3;
    double t345 = t22*t23*t45*3.064210757541298E-3;
    double t346 = t9*t10*t35*-3.104080344633556E-3;
    double t347 = t9*t10*t35*3.104080344633556E-3;
    double t348 = t9*t11*t34*3.104080344633556E-3;
    double t349 = t10*t11*t33*3.104080344633556E-3;
    double t350 = t21*t22*t47*3.104080344633556E-3;
    double t351 = t21*t23*t46*3.104080344633556E-3;
    double t352 = t22*t23*t45*3.104080344633556E-3;
    double t353 = t6*t7*t32*3.105990081579729E-3;
    double t354 = t6*t8*t31*3.105990081579729E-3;
    double t355 = t7*t8*t30*3.105990081579729E-3;
    double t356 = t18*t19*t44*-3.105990081579729E-3;
    double t357 = t18*t19*t44*3.105990081579729E-3;
    double t358 = t18*t20*t43*-3.105990081579729E-3;
    double t359 = t18*t20*t43*3.105990081579729E-3;
    double t360 = t19*t20*t42*3.105990081579729E-3;
    double t361 = t6*t31*t32*3.397508858570615E-1;
    double t362 = t7*t30*t32*3.397508858570615E-1;
    double t363 = t8*t30*t31*3.397508858570615E-1;
    double t364 = t18*t43*t44*3.397508858570615E-1;
    double t365 = t19*t42*t44*3.397508858570615E-1;
    double t366 = t20*t42*t43*-3.397508858570615E-1;
    double t367 = t20*t42*t43*3.397508858570615E-1;
    double t368 = t6*t31*t32*3.399783924207052E-1;
    double t369 = t7*t30*t32*3.399783924207052E-1;
    double t370 = t8*t30*t31*3.399783924207052E-1;
    double t371 = t18*t43*t44*3.399783924207052E-1;
    double t372 = t19*t42*t44*3.399783924207052E-1;
    double t373 = t20*t42*t43*-3.399783924207052E-1;
    double t374 = t20*t42*t43*3.399783924207052E-1;
    double t375 = t33*t34*t35*-2.87566285868891E-1;
    double t376 = t33*t34*t35*2.87566285868891E-1;
    double t377 = t45*t46*t47*-2.87566285868891E-1;
    double t378 = t45*t46*t47*2.87566285868891E-1;
    double t379 = t33*t34*t35*-2.879825211612492E-1;
    double t380 = t33*t34*t35*2.879825211612492E-1;
    double t381 = t45*t46*t47*2.879825211612492E-1;
    double t382 = t6*t31*t32*3.020283789547073E-3;
    double t383 = t7*t30*t32*3.020283789547073E-3;
    double t384 = t8*t30*t31*-3.020283789547073E-3;
    double t385 = t8*t30*t31*3.020283789547073E-3;
    double t386 = t18*t43*t44*3.020283789547073E-3;
    double t387 = t19*t42*t44*3.020283789547073E-3;
    double t388 = t20*t42*t43*3.020283789547073E-3;
    double t389 = t9*t34*t35*3.064210757541298E-3;
    double t390 = t10*t33*t35*3.064210757541298E-3;
    double t391 = t11*t33*t34*-3.064210757541298E-3;
    double t392 = t11*t33*t34*3.064210757541298E-3;
    double t393 = t21*t46*t47*-3.064210757541298E-3;
    double t394 = t21*t46*t47*3.064210757541298E-3;
    double t395 = t22*t45*t47*3.064210757541298E-3;
    double t396 = t23*t45*t46*-3.064210757541298E-3;
    double t397 = t23*t45*t46*3.064210757541298E-3;
    double t398 = t9*t34*t35*3.104080344633556E-3;
    double t399 = t10*t33*t35*3.104080344633556E-3;
    double t400 = t11*t33*t34*3.104080344633556E-3;
    double t401 = t21*t46*t47*3.104080344633556E-3;
    double t402 = t22*t45*t47*-3.104080344633556E-3;
    double t403 = t22*t45*t47*3.104080344633556E-3;
    double t404 = t23*t45*t46*-3.104080344633556E-3;
    double t405 = t23*t45*t46*3.104080344633556E-3;
    double t406 = t6*t31*t32*3.105990081579729E-3;
    double t407 = t7*t30*t32*3.105990081579729E-3;
    double t408 = t8*t30*t31*-3.105990081579729E-3;
    double t409 = t8*t30*t31*3.105990081579729E-3;
    double t410 = t18*t43*t44*3.105990081579729E-3;
    double t411 = t19*t42*t44*3.105990081579729E-3;
    double t412 = t20*t42*t43*3.105990081579729E-3;
    double t413 = t30*t31*t32*3.397508858570615E-1;
    double t414 = t42*t43*t44*-3.397508858570615E-1;
    double t415 = t42*t43*t44*3.397508858570615E-1;
    double t416 = t30*t31*t32*3.399783924207052E-1;
    double t417 = t42*t43*t44*-3.399783924207052E-1;
    double t418 = t42*t43*t44*3.399783924207052E-1;
    double t419 = t30*t31*t32*-3.020283789547073E-3;
    double t420 = t30*t31*t32*3.020283789547073E-3;
    double t421 = t42*t43*t44*3.020283789547073E-3;
    double t422 = t33*t34*t35*3.064210757541298E-3;
    double t423 = t45*t46*t47*-3.064210757541298E-3;
    double t424 = t45*t46*t47*3.064210757541298E-3;
    double t425 = t33*t34*t35*3.104080344633556E-3;
    double t426 = t45*t46*t47*-3.104080344633556E-3;
    double t427 = t45*t46*t47*3.104080344633556E-3;
    double t428 = t30*t31*t32*-3.105990081579729E-3;
    double t429 = t30*t31*t32*3.105990081579729E-3;
    double t430 = t42*t43*t44*3.105990081579729E-3;
    double t431 = t70+t101;
    double t432 = t71+t103;
    double t433 = t100+t151;
    double t434 = t102+t153;
    double t435 = t126+t154;
    double t436 = t129+t155;
    double t437 = t64+t239;
    double t438 = t65+t241;
    double t439 = t84+t212;
    double t440 = t87+t213;
    double t441 = t125+t225;
    double t442 = t128+t226;
    double t443 = t143+t217;
    double t444 = t146+t218;
    double t445 = t137+t238;
    double t446 = t139+t240;
    double t447 = t83+t251;
    double t448 = t86+t252;
    double t449 = t142+t256;
    double t450 = t145+t257;
    double t486 = t66+t111+t170+t220;
    double t487 = t76+t88+t185+t201;
    double t488 = t90+t121+t162+t183;
    double t489 = t112+t130+t140+t168;
    double t490 = t67+t116+t172+t223;
    double t491 = t78+t92+t191+t203;
    double t492 = t93+t122+t165+t187;
    double t493 = t120+t132+t141+t171;
    double t494 = t68+t109+t177+t219;
    double t495 = t75+t94+t186+t210;
    double t496 = t95+t123+t163+t182;
    double t497 = t113+t131+t147+t176;
    double t498 = t69+t117+t180+t222;
    double t499 = t77+t97+t190+t211;
    double t500 = t99+t124+t166+t189;
    double t501 = t119+t133+t149+t178;
    double t502 = t61+t134+t231+t232;
    double t503 = t79+t80+t198+t248;
    double t504 = t63+t135+t233+t234;
    double t505 = t81+t82+t199+t250;
    double t547 = t272+t379+t390+t398;
    double t548 = t280+t375+t389+t399;
    double t549 = t318+t324+t339+t425;
    double t550 = t317+t325+t346+t422;
    double t551 = t276+t381+t395+t401;
    double t552 = t284+t378+t394+t403;
    double t553 = t322+t328+t343+t426;
    double t554 = t321+t329+t350+t423;
    double t556 = t301+t382+t407+t416;
    double t557 = t309+t383+t406+t413;
    double t558 = t331+t361+t369+t428;
    double t559 = t353+t362+t368+t419;
    double t560 = t306+t386+t411+t417;
    double t561 = t314+t387+t410+t414;
    double t562 = t334+t364+t372+t430;
    double t563 = t356+t365+t371+t421;
    double t565 = t274+t283+t295+t392;
    double t566 = t275+t281+t292+t400;
    double t567 = t262+t319+t341+t349;
    double t568 = t259+t326+t342+t348;
    double t569 = t278+t287+t297+t396;
    double t570 = t279+t286+t294+t404;
    double t571 = t263+t323+t344+t352;
    double t572 = t260+t330+t345+t351;
    double t574 = t289+t304+t313+t408;
    double t575 = t298+t305+t311+t384;
    double t576 = t265+t332+t355+t370;
    double t577 = t268+t333+t354+t363;
    double t578 = t290+t307+t316+t412;
    double t579 = t299+t308+t315+t388;
    double t580 = t267+t336+t360+t373;
    double t581 = t271+t338+t359+t366;
    double t662 = t235+t243+t273+t282+t296+t391;
    double t667 = t175+t247+t285+t377+t393+t402;
    double t670 = t159+t194+t288+t303+t312+t409;
    double t671 = t107+t230+t270+t337+t358+t367;
    double t451 = t431*x10;
    double t452 = t431*x11;
    double t453 = t432*x25;
    double t454 = t432*x26;
    double t455 = t433*x10;
    double t456 = t434*x25;
    double t457 = t435*x15;
    double t458 = t435*x16;
    double t459 = t436*x30;
    double t460 = t436*x31;
    double t461 = t437*x10;
    double t462 = t437*x11;
    double t463 = t438*x25;
    double t464 = t438*x26;
    double t465 = t439*x10;
    double t466 = t439*x11;
    double t467 = t440*x25;
    double t468 = t440*x26;
    double t469 = t441*x15;
    double t470 = t442*x30;
    double t471 = t443*x18;
    double t472 = t443*x19;
    double t473 = t444*x33;
    double t474 = t444*x34;
    double t475 = t445*x10;
    double t476 = t446*x25;
    double t477 = t447*x10;
    double t478 = t448*x25;
    double t479 = t449*x18;
    double t480 = t450*x33;
    double t506 = t486*x20;
    double t507 = t488*x20;
    double t508 = t486*x21;
    double t509 = t487*x21;
    double t510 = t488*x21;
    double t511 = t489*x21;
    double t512 = t490*x35;
    double t513 = t492*x35;
    double t514 = t490*x36;
    double t515 = t491*x36;
    double t516 = t492*x36;
    double t517 = t493*x36;
    double t518 = t494*x20;
    double t519 = t496*x20;
    double t520 = t494*x21;
    double t521 = t495*x21;
    double t522 = t496*x21;
    double t523 = t497*x21;
    double t524 = t498*x35;
    double t525 = t500*x35;
    double t526 = t498*x36;
    double t527 = t499*x36;
    double t528 = t500*x36;
    double t529 = t501*x36;
    double t530 = t502*x12;
    double t531 = t503*x12;
    double t532 = t502*x13;
    double t533 = t503*x13;
    double t534 = t504*x27;
    double t535 = t505*x27;
    double t536 = t504*x28;
    double t537 = t505*x28;
    double t584 = t566*x17;
    double t585 = t568*x17;
    double t586 = t565*x18;
    double t587 = t566*x18;
    double t588 = t567*x18;
    double t589 = t568*x18;
    double t590 = t566*x19;
    double t591 = t568*x19;
    double t592 = t570*x32;
    double t593 = t572*x32;
    double t594 = t569*x33;
    double t595 = t570*x33;
    double t596 = t571*x33;
    double t597 = t572*x33;
    double t598 = t570*x34;
    double t599 = t572*x34;
    double t600 = t575*x14;
    double t601 = t576*x14;
    double t602 = t574*x15;
    double t603 = t575*x15;
    double t604 = t576*x15;
    double t605 = t577*x15;
    double t606 = t575*x16;
    double t607 = t576*x16;
    double t608 = t579*x29;
    double t609 = t580*x29;
    double t610 = t578*x30;
    double t611 = t579*x30;
    double t612 = t580*x30;
    double t613 = t581*x30;
    double t614 = t579*x31;
    double t615 = t580*x31;
    double t616 = t547*x17;
    double t617 = t549*x17;
    double t618 = t547*x18;
    double t619 = t548*x18;
    double t620 = t549*x18;
    double t621 = t550*x18;
    double t622 = t547*x19;
    double t623 = t549*x19;
    double t624 = t551*x32;
    double t625 = t553*x32;
    double t626 = t551*x33;
    double t627 = t552*x33;
    double t628 = t553*x33;
    double t629 = t554*x33;
    double t630 = t551*x34;
    double t631 = t553*x34;
    double t632 = t556*x14;
    double t633 = t559*x14;
    double t634 = t556*x15;
    double t635 = t557*x15;
    double t636 = t558*x15;
    double t637 = t559*x15;
    double t638 = t556*x16;
    double t639 = t559*x16;
    double t640 = t560*x29;
    double t641 = t563*x29;
    double t642 = t560*x30;
    double t643 = t561*x30;
    double t644 = t562*x30;
    double t645 = t563*x30;
    double t646 = t560*x31;
    double t647 = t563*x31;
    double t657 = t72+t195+t557;
    double t658 = t106+t156+t558;
    double t659 = t74+t197+t561;
    double t660 = t108+t157+t562;
    double t661 = t205+t254+t567;
    double t663 = t208+t255+t571;
    double t664 = t236+t246+t569;
    double t665 = t174+t244+t548;
    double t666 = t207+t214+t550;
    double t668 = t209+t215+t554;
    double t669 = t105+t228+t577;
    double t672 = t160+t196+t578;
    double t677 = t667*x32;
    double t679 = t667*x34;
    double t682 = t670*x16;
    double t683 = t671*x31;
    double t694 = t662*x19;
    double t481 = -t456;
    double t482 = -t469;
    double t483 = -t476;
    double t484 = -t477;
    double t485 = -t479;
    double t538 = -t511;
    double t539 = -t517;
    double t540 = -t521;
    double t541 = -t527;
    double t542 = t474+t480;
    double t543 = t452+t455;
    double t544 = t460+t470;
    double t545 = t462+t475;
    double t546 = t468+t478;
    double t648 = -t591;
    double t649 = -t598;
    double t650 = -t606;
    double t651 = -t615;
    double t652 = -t623;
    double t653 = -t626;
    double t654 = -t630;
    double t655 = -t638;
    double t656 = -t647;
    double t673 = t665*x17;
    double t674 = t666*x17;
    double t675 = t665*x19;
    double t676 = t666*x19;
    double t678 = t668*x32;
    double t680 = t668*x34;
    double t681 = t669*x16;
    double t684 = t672*x31;
    double t685 = t657*x14;
    double t686 = t658*x14;
    double t687 = t657*x16;
    double t688 = t658*x16;
    double t689 = t659*x29;
    double t690 = t660*x29;
    double t691 = t659*x31;
    double t692 = t660*x31;
    double t693 = t661*x19;
    double t695 = t663*x34;
    double t696 = t664*x34;
    double t699 = -t694;
    double t700 = t507+t509;
    double t701 = t513+t515;
    double t702 = t518+t523;
    double t703 = t524+t529;
    double t704 = t530+t532;
    double t705 = t531+t533;
    double t706 = t534+t536;
    double t707 = t535+t537;
    double t717 = t584+t586+t622;
    double t718 = t593+t596+t631;
    double t719 = t590+t616+t619;
    double t720 = t599+t625+t629;
    double t721 = t601+t605+t639;
    double t722 = t608+t610+t646;
    double t723 = t607+t633+t636;
    double t724 = t614+t640+t643;
    double t555 = t472+t485;
    double t564 = t454+t481;
    double t573 = t458+t482;
    double t582 = t464+t483;
    double t583 = t466+t484;
    double t697 = -t676;
    double t698 = -t693;
    double t708 = t506+t538;
    double t709 = t512+t539;
    double t710 = t519+t540;
    double t711 = t525+t541;
    double t712 = t704*y13;
    double t713 = t705*y13;
    double t714 = t706*y28;
    double t715 = t707*y28;
    double t725 = t585+t588+t652;
    double t726 = t592+t594+t654;
    double t727 = t617+t621+t648;
    double t728 = t624+t627+t649;
    double t729 = t600+t602+t655;
    double t730 = t609+t613+t656;
    double t731 = t632+t635+t650;
    double t732 = t641+t644+t651;
    double t733 = t628+t678+t695;
    double t734 = t637+t681+t686;
    double t735 = t634+t682+t685;
    double t736 = t645+t683+t690;
    double t737 = t642+t684+t689;
    double t738 = t618+t673+t699;
    double t740 = t653+t677+t696;
    double t716 = -t713;
    double t739 = t620+t674+t698;

    Jxy_partial_dq.setZero();
    Jxy_partial_dq(0, 28) = y29*(t611+t691+x29*(t18*-5.699060997402858E-2+t42*1.034589188110661E-3+t672))+t722*y30+t737*y31;
    Jxy_partial_dq(0, 29) = y30*(t611+t561*x31+t578*x29)+t722*y29+t724*y31;
    Jxy_partial_dq(0, 30) = t724*y30+t737*y29+y31*(t611+t691+t672*x29);
    Jxy_partial_dq(0, 34) = y35*(t514-x35*(t55+t58+t493))+t709*y36;
    Jxy_partial_dq(0, 35) = t709*y35+y36*(t514-t493*x35);
    Jxy_partial_dq(1, 28) = y29*(t612-t692+x29*(t18*1.034589188110661E-3+t42*5.699060997402858E-2-t107-t230+t581))+t730*y30-t736*y31;
    Jxy_partial_dq(1, 29) = y30*(t612-t562*x31+t581*x29)+t730*y29-t732*y31;
    Jxy_partial_dq(1, 30) = -t732*y30-t736*y29-y31*(-t612+t692+t671*x29);
    Jxy_partial_dq(1, 34) = t701*y36+y35*(t516-x35*(t51+t60-t78-t92+t192+t204));
    Jxy_partial_dq(1, 35) = t701*y35+y36*(t516+t491*x35);
    Jxy_partial_dq(2, 29) = t544*y31+y30*(t459+t442*x31);
    Jxy_partial_dq(2, 30) = t544*y30+y31*(t459-x31*(t20*3.39756885202024E-1-t128+t227));
    Jxy_partial_dq(2, 35) = x36*y36*(t49*1.79E-2+t53);
    Jxy_partial_dq(3, 31) = y32*(t595+t679+x32*(t21*5.699060997402856E-2+t45*1.034589188110661E-3+t664))+t726*y33+t740*y34;
    Jxy_partial_dq(3, 32) = y33*(t595-t552*x34+t569*x32)+t726*y32-t728*y34;
    Jxy_partial_dq(3, 33) = -t728*y33+t740*y32+y34*(t595+t679+t664*x32);
    Jxy_partial_dq(3, 34) = -y35*(t526+x35*(t55+t58+t501))-t703*y36;
    Jxy_partial_dq(3, 35) = -t703*y35-y36*(t526+t501*x35);
    Jxy_partial_dq(4, 31) = y32*(t597+t680+x32*(t21*-1.034589188110661E-3+t45*5.699060997402856E-2+t663))+t718*y33+t733*y34;
    Jxy_partial_dq(4, 32) = y33*(t597+t554*x34+t571*x32)+t718*y32+t720*y34;
    Jxy_partial_dq(4, 33) = t720*y33+t733*y32+y34*(t597+t680+t663*x32);
    Jxy_partial_dq(4, 34) = -t711*y36-y35*(t528+x35*(t51+t60-t77+t98+t193-t211));
    Jxy_partial_dq(4, 35) = -t711*y35-y36*(t528-t499*x35);
    Jxy_partial_dq(5, 32) = -t542*y34-y33*(t473+t450*x34);
    Jxy_partial_dq(5, 33) = -y34*(t473-x34*(t23*2.875818595898751E-1-t145+t258))-t542*y33;
    Jxy_partial_dq(5, 35) = -x36*y36*(t49*1.81E-2-t53);
    Jxy_partial_dq(6, 13) = y14*(t603-t687+x14*(t6*5.699060997402858E-2+t30*1.034589188110661E-3-t159-t194+t574))+t729*y15-t735*y16;
    Jxy_partial_dq(6, 14) = y15*(t603-t557*x16+t574*x14)+t729*y14-t731*y16;
    Jxy_partial_dq(6, 15) = -t731*y15-t735*y14-y16*(-t603+t687+t670*x14);
    Jxy_partial_dq(6, 19) = -y20*(t508-x20*(t54+t57+t489))-t708*y21;
    Jxy_partial_dq(6, 20) = -t708*y20-y21*(t508-t489*x20);
    Jxy_partial_dq(7, 13) = -y14*(t604+t688+x14*(t6*-1.034589188110661E-3+t30*5.699060997402858E-2+t669))-t721*y15-t734*y16;
    Jxy_partial_dq(7, 14) = -y15*(t604+t558*x16+t577*x14)-t721*y14-t723*y16;
    Jxy_partial_dq(7, 15) = -t723*y15-t734*y14-y16*(t604+t688+t669*x14);
    Jxy_partial_dq(7, 19) = t700*y21+y20*(t510-x20*(t50-t59-t76+t89-t185+t202));
    Jxy_partial_dq(7, 20) = t700*y20+y21*(t510+t487*x20);
    Jxy_partial_dq(8, 14) = -t573*y16-y15*(t457-t441*x16);
    Jxy_partial_dq(8, 15) = -t573*y15-y16*(t457+x16*(t8*3.39756885202024E-1-t441));
    Jxy_partial_dq(8, 20) = x21*y21*(t37*1.79E-2-t52);
    Jxy_partial_dq(9, 16) = y17*(t587+t675-x17*(t9*5.699060997402856E-2-t33*1.034589188110661E-3+t662))+t717*y18+t738*y19;
    Jxy_partial_dq(9, 17) = y18*(t587+t548*x19+t565*x17)+t717*y17+t719*y19;
    Jxy_partial_dq(9, 18) = t719*y18+t738*y17+y19*(t587+t675-t662*x17);
    Jxy_partial_dq(9, 19) = y20*(t520+x20*(t54+t57+t497))+t702*y21;
    Jxy_partial_dq(9, 20) = t702*y20+y21*(t520+t497*x20);
    Jxy_partial_dq(10, 16) = -y17*(t589+t697+x17*(t9*1.034589188110661E-3+t33*5.699060997402856E-2+t661))-t725*y18+t739*y19;
    Jxy_partial_dq(10, 17) = -y18*(t589-t550*x19+t567*x17)-t725*y17+t727*y19;
    Jxy_partial_dq(10, 18) = t727*y18+t739*y17-y19*(t589+t697+t661*x17);
    Jxy_partial_dq(10, 19) = -t710*y21-y20*(t522-x20*(-t50+t59+t495));
    Jxy_partial_dq(10, 20) = -t710*y20-y21*(t522-t495*x20);
    Jxy_partial_dq(11, 17) = t555*y19+y18*(t471-t449*x19);
    Jxy_partial_dq(11, 18) = t555*y18+y19*(t471+x19*(t11*2.875818595898751E-1-t449));
    Jxy_partial_dq(11, 20) = -x21*y21*(t37*1.81E-2+t52);
    Jxy_partial_dq(12, 24) = t546*y26+y25*(t467+t448*x26);
    Jxy_partial_dq(12, 25) = t546*y25+y26*(t467+x26*(t15*4.365860704565494E-4+t448));
    Jxy_partial_dq(12, 26) = t714+y27*(t536-x27*(t16*4.95436E-1+t40*6.740600000000002E-2-t63+t136-t233-t234));
    Jxy_partial_dq(12, 27) = t714+t706*y27;
    Jxy_partial_dq(13, 24) = -t582*y26-y25*(t463-t446*x26);
    Jxy_partial_dq(13, 25) = -t582*y25-y26*(t463-x26*(t15*3.566153386244494E-2+t446));
    Jxy_partial_dq(13, 26) = t715+y27*(t537+x27*(t16*6.740600000000002E-2-t40*4.95436E-1+t505));
    Jxy_partial_dq(13, 27) = t715+t707*y27;
    Jxy_partial_dq(14, 24) = t564*y26+y25*(t453-t434*x26);
    Jxy_partial_dq(14, 25) = y26*(t453+x26*(t15*4.987264424463382E-1-t434))+t564*y25;
    Jxy_partial_dq(15, 9) = -t583*y11-y10*(t465-t447*x11);
    Jxy_partial_dq(15, 10) = -t583*y10-y11*(t465-x11*(t3*4.365860704565494E-4+t447));
    Jxy_partial_dq(15, 11) = t712+y12*(t532+x12*(t4*4.95436E-1-t28*6.740600000000002E-2+t502));
    Jxy_partial_dq(15, 12) = t712+t704*y12;
    Jxy_partial_dq(16, 9) = -t545*y11-y10*(t461+t445*x11);
    Jxy_partial_dq(16, 10) = -t545*y10-y11*(t461+x11*(t3*3.566153386244494E-2+t445));
    Jxy_partial_dq(16, 11) = t716-y12*(t533-x12*(t4*6.740600000000002E-2+t28*4.95436E-1-t79-t80-t198+t249));
    Jxy_partial_dq(16, 12) = t716-t705*y12;
    Jxy_partial_dq(17, 9) = -t543*y11-y10*(t451+t433*x11);
    Jxy_partial_dq(17, 10) = -t543*y10-y11*(t451-x11*(t3*4.987264424463382E-1-t100+t152));

    fkPtr_->compute(0, 
                    contact_joint_id, 
                    q, 
                    nullptr, 
                    &stance_foot_endT, 
                    3);

    fkPtr_->getTranslationThirdOrderTensor(x, TOx_translation);
    fkPtr_->getRPYThirdOrderTensor(x, TOx_rotation);

    Jxy_partial_dq.row(18) = TOx_translation(0) * y;
    Jxy_partial_dq.row(19) = TOx_translation(1) * y;
    Jxy_partial_dq.row(20) = TOx_translation(2) * y;
    Jxy_partial_dq.row(21) = TOx_rotation(0) * y;
    Jxy_partial_dq.row(22) = TOx_rotation(1) * y;
    Jxy_partial_dq.row(23) = TOx_rotation(2) * y;
}

}; // namespace Digit
}; // namespace RAPTOR
