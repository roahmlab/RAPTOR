#include "DigitConstrainedInverseDynamics.h"

namespace IDTO {
namespace Digit {

DigitConstrainedInverseDynamics::DigitConstrainedInverseDynamics(const Model& model_input, 
                                                                 int N_input, 
                                                                 int numDependentJoints_input,
                                                                 std::unique_ptr<DynamicsConstraints>& dynamics_constraints_input) :
    ConstrainedInverseDynamics(model_input, N_input, numDependentJoints_input, dynamics_constraints_input) {

}

void DigitConstrainedInverseDynamics::setupJointPosition(VecX& q) {
    VecX qcopy = q;

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
    qcopy(modelPtr_->getJointId("right_toe_A_rod") - 1) = -asin(t2*5.114967905744713E-1+t3*4.072902857734452E-1+t4*5.15019991000208E-1-t5*7.189052210514619E-1-t8*1.914308797309167E-1+t9*1.020482430826409E-1+t10*2.762037569148706E-1+t11*4.461078267300372E-1-t2*t3*5.669643819573642E-1+t2*t5*9.826313004441752E-1+t3*t4*7.066611122892024E-1-t4*t5*8.751268427354268E-2-t2*t9*1.526013089230614E-1+t3*t8*2.146475284950957E-1-t2*t11*6.373722738059316E-1-t3*t10*4.061233553924697E-1-t4*t9*2.591685269591633E-1-t5*t8*3.141536266927853E-1+t4*t11*1.154141052196728E-1+t5*t10*5.111686735370716E-2+t8*t9*2.629326893467769E-2+t8*t11*2.261999172185061E-1+t9*t10*1.714455893661597E-1-t10*t11*5.322298539488757E-2-3.496118766033006E-1);
    qcopy(modelPtr_->getJointId("right_A2") - 1) = asin(t2*8.08508988653757+t3*8.033204820099098+t4*1.820663506981432+t5*9.760946687487468E-1-t8*2.143120822774579-t9*1.998293355526279-t10*1.095961671709748-t11*4.953901287885355E-1-t2*t3*1.273036938982479E+1-t2*t5*1.117825053271672-t3*t4*2.420405429552186+t4*t5*1.170037217867019+t2*t9*2.835751824751641+t3*t8*3.021210259875125+t2*t11*5.805706025081291E-1+t3*t10*1.468854423683356+t4*t9*5.420356015832145E-1+t5*t8*2.326217949887465E-1+t4*t11*1.600391744429722E-1+t5*t10*1.170908991628156E-1-t8*t9*4.40337591097512E-1-t8*t11*1.065731420572547E-1-t9*t10*3.400548068856727E-1-t10*t11*4.121086576936894E-1-4.663948776295725);
    qcopy(modelPtr_->getJointId("right_toe_B_rod") - 1) = -asin(t2*(-6.277469604269263E-1)-t3*6.676022839778379E-1-t4*8.399423316628958E-1+t5*1.529719504893023E-1+t8*1.44381731007234E-1+t9*4.870768316848186E-1+t10*5.196004288358185E-1+t11*4.899787657102767E-1+t2*t3*1.038356653790214+t2*t5*1.224128324227442+t3*t4*1.136249345040226-t4*t5*2.962244825741616E-1-t2*t9*7.061936994304631E-1-t3*t8*2.508516890107657E-1-t2*t11*7.147313228957868E-1-t3*t10*7.361945978788478E-1-t4*t9*3.638713473685915E-1-t5*t8*4.137119943797418E-1+t4*t11*2.198247731883075E-1+t5*t10*1.190720324216782E-1+t8*t9*1.91374695742122E-1+t8*t11*2.674233531253601E-1+t9*t10*2.592073581826719E-1-t10*t11*1.033333326963321E-1+3.90361533934657E-1);
    qcopy(modelPtr_->getJointId("right_B2") - 1) = asin(t2*9.834460864042423+t3*1.012210878062529E+1+t4*4.114199768067971+t5*2.918615505121678-t8*2.33635398066763-t9*2.824506662141518-t10*2.436915399445716-t11*1.657161229950764-t2*t3*1.56824564483544E+1-t2*t5*4.071020737587724-t3*t4*5.741194151827515+t4*t5*1.412990674523749+t2*t9*3.933187185458805+t3*t8*3.362495208957303+t2*t11*2.232327688060054+t3*t10*3.377470435585635+t4*t9*1.500870022094839+t5*t8*1.104866050029796+t4*t11*1.111537011807542E-1+t5*t10*1.885808892070127E-1-t8*t9*5.75161230009884E-1-t8*t11*6.064451339600673E-1-t9*t10*9.059551372991537E-1-t10*t11*4.795029892364048E-1-5.829415122169579);
    qcopy(modelPtr_->getJointId("right_toe_pitch") - 1) = -asin(t2*1.893989559445235E+1+t3*1.920869550741533E+1+t4*6.868951098405778+t5*3.51561215609706-t8*4.080342460091422-t9*4.456343821978674-t10*3.407846365996183-t11*2.503419572692307-t2*t3*2.929974638198693E+1-t2*t5*5.157359318839895-t3*t4*8.730366732237334+t4*t5*1.967066798625298+t2*t9*5.895002722892592+t3*t8*5.481616764469075+t2*t11*3.209642327668312+t3*t10*4.724342470672166+t4*t9*2.235861962074205+t5*t8*1.242934577497499+t4*t11*1.119458771782081+t5*t10*1.139096075125015-t8*t9*3.306304563913454E-1-t8*t11*7.650487856200098E-1-t9*t10*1.240832663118515-t10*t11*1.549291982370939-1.135688467845833E+1);
    qcopy(modelPtr_->getJointId("right_toe_roll") - 1) = -asin(t2*1.925184467728597+t3*7.345366777287639+t4*2.843699508105022E+1+t5*2.720900626308319E+1+t8*3.57433462778918-t9*6.332457856504018-t10*1.691966583575934E+1-t11*1.60681725794845E+1-t2*t3*6.299462363896048-t2*t5*4.117199083152367E+1-t3*t4*4.282360516873334E+1+t4*t5*7.910198316638063E-2+t2*t9*7.854243405249483-t3*t8*4.011841131466795+t2*t11*2.320584908063003E+1+t3*t10*2.435545111239733E+1+t4*t9*1.284334673751474E+1+t5*t8*1.242519321726289E+1-t4*t11*1.106503053903957+t5*t10*8.20592115289103E-1-t8*t9*6.089534825702931E-1-t8*t11*7.206363286390523-t9*t10*7.502224690460121+t10*t11*1.255468934813467E-1-3.446675672448973);

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
    qcopy(modelPtr_->getJointId("left_toe_A_rod") - 1) = -asin(t2*(-5.114926838701811E-1)-t3*4.072861778028993E-1+t4*5.150198314974532E-1-t5*7.189053388109297E-1+t8*1.914297694169048E-1-t9*1.020493568564005E-1+t10*2.762038362432678E-1+t11*4.461078824742524E-1+t2*t3*5.669587551620541E-1+t2*t5*9.826314482751646E-1+t3*t4*7.066613185154341E-1+t4*t5*8.751254943270917E-2+t2*t9*1.526028367712922E-1-t3*t8*2.146460048968412E-1-t2*t11*6.3737234338624E-1-t3*t10*4.061234579946822E-1-t4*t9*2.591685735898548E-1-t5*t8*3.141536565807085E-1-t4*t11*1.15414032157577E-1-t5*t10*5.111679397817798E-2-t8*t9*2.629368451925844E-2+t8*t11*2.261999309581469E-1+t9*t10*1.714456126031132E-1+t10*t11*5.322294548977349E-2+3.49608876912165E-1);
    qcopy(modelPtr_->getJointId("left_A2") - 1) = asin(t2*8.085018160207163+t3*8.033133064086231-t4*1.820662389025083-t5*9.760941431166779E-1-t8*2.143101587807895-t9*1.998274051537184+t10*1.095961134647075+t11*4.953899201857669E-1-t2*t3*1.273027136495291E+1+t2*t5*1.11782444127773+t3*t4*2.420403987960523+t4*t5*1.17003991016308+t2*t9*2.835725404415498+t3*t8*3.021183922843985-t2*t11*5.805703711334091E-1-t3*t10*1.468853731639082-t4*t9*5.420352775054091E-1-t5*t8*2.326217084069311E-1+t4*t11*1.600377315075955E-1+t5*t10*1.170894522681423E-1-t8*t9*4.403304524513951E-1+t8*t11*1.065731191784157E-1+t9*t10*3.400546515442943E-1-t10*t11*4.121078791939159E-1-4.663896238435751);
    qcopy(modelPtr_->getJointId("left_toe_B_rod") - 1) = -asin(t2*6.277394162562625E-1+t3*6.6759473381467E-1-t4*8.399421448959555E-1+t5*1.52972128181329E-1-t8*1.443797201746081E-1-t9*4.870748185116889E-1+t10*5.196003431996498E-1+t11*4.899786813403554E-1-t2*t3*1.038346357873659+t2*t5*1.224128102824049+t3*t4*1.13624910964581+t4*t5*2.962247949009167E-1+t2*t9*7.06190949266457E-1+t3*t8*2.508489399062072E-1-t2*t11*7.147312184882999E-1-t3*t10*7.361944907642444E-1-t4*t9*3.638712999144287E-1-t5*t8*4.137119506580644E-1-t4*t11*2.198249429188003E-1-t5*t10*1.190722007592886E-1-t8*t9*1.913739574570494E-1+t8*t11*2.67423333074767E-1+t9*t10*2.592073373292207E-1+t10*t11*1.033334244581574E-1-3.903559984881447E-1);
    qcopy(modelPtr_->getJointId("left_B2") - 1) = asin(t2*9.83436298937519+t3*1.012201086685809E+1-t4*4.114198671144258-t5*2.918615149671792-t8*2.3363277301986-t9*2.824480323210008+t10*2.436914889769644+t11*1.657161129062592-t2*t3*1.568232268107518E+1+t2*t5*4.071020370775729+t3*t4*5.741192746046224+t4*t5*1.412994342972937+t2*t9*3.933151135492005+t3*t8*3.362459265404445-t2*t11*2.232327610328257-t3*t10*3.377469784128673-t4*t9*1.500869712732554-t5*t8*1.104866037480153+t4*t11*1.111517321363341E-1+t5*t10*1.885789156342841E-1-t8*t9*5.751514893646594E-1+t8*t11*6.064451565041576E-1+t9*t10*9.05954995114355E-1-t10*t11*4.795019258370027E-1-5.829343436697216);
    qcopy(modelPtr_->getJointId("left_toe_pitch") - 1) = asin(t2*1.893971200186481E+1+t3*1.920851184673561E+1-t4*6.868948549897667-t5*3.515611050123151-t8*4.080293186785779-t9*4.456294384490135+t10*3.407845158508307+t11*2.503419165594914-t2*t3*2.92994954151725E+1+t2*t5*5.157358065307949+t3*t4*8.730363458834317+t4*t5*1.967073609987012+t2*t9*5.894935044834759+t3*t8*5.481549285471837-t2*t11*3.209641901519396-t3*t10*4.72434092223045-t4*t9*2.235861236273534-t5*t8*1.242934428625482+t4*t11*1.119455114838831+t5*t10*1.139092408157223-t8*t9*3.306121591996909E-1+t8*t11*7.650487658943779E-1+t9*t10*1.240832321393588-t10*t11*1.549290006003741-1.135675024134487E+1);
    qcopy(modelPtr_->getJointId("left_toe_roll") - 1) = asin(t2*1.92500223858939+t3*7.345184505784512-t4*2.843699897601034E+1-t5*2.720901076343499E+1+t8*3.574383636815427-t9*6.332408732694388+t10*1.691966794969343E+1+t11*1.606817500398276E+1-t2*t3*6.299213096460508+t2*t5*4.117199681273369E+1+t3*t4*4.282361029947926E+1+t4*t5*7.910853681196074E-2+t2*t9*7.854176116345911-t3*t8*4.011908279005063-t2*t11*2.320585230931768E+1-t3*t10*2.435545390288711E+1-t4*t9*1.284334796900915E+1-t5*t8*1.242519468406818E+1-t4*t11*1.106506602910529+t5*t10*8.205885613072367E-1-t8*t9*6.089352667535674E-1+t8*t11*7.206364083227697+t9*t10*7.502225364895965+t10*t11*1.255488247882333E-1-3.446542349892159);

    // right knee close loop
    q1 = q(modelPtr_->getJointId("right_knee") - 1);
    t2 = cos(q1);
    t3 = sin(q1);
    t4 = q1*2.0;
    t5 = cos(t4);
    t6 = sin(t4);
    qcopy(modelPtr_->getJointId("right_tarsus") - 1) = -asin(t2*1.155848969647063E-3+t3*1.004686948291003+t5*1.274417498011625E-4-t6*1.785981355062532E-3-1.132590494159057E-2);
    qcopy(modelPtr_->getJointId("right_achilles_rod") - 1) = -asin(t2*(-1.587289102030986E-3)-t3*1.001736520672665+t5*3.407131509821247E-4+t6*9.646678263881318E-4+1.539911054998293E-3);
    qcopy(modelPtr_->getJointId("right_ach2") - 1) = asin(t2*(-7.197863326636346E-2)-t3*8.929579539511397E-3+t5*2.669904889172627E-4+t6*8.46571305589265E-5+7.18964949007849E-2);

    // left knee close loop
    q1 = q(modelPtr_->getJointId("left_knee") - 1);
    t2 = cos(q1);
    t3 = sin(q1);
    t4 = q1*2.0;
    t5 = cos(t4);
    t6 = sin(t4);
    qcopy(modelPtr_->getJointId("left_tarsus") - 1) = asin(t2*1.155848972188414E-3-t3*1.004686948291033+t5*1.274417489907877E-4+t6*1.785981355072367E-3-1.132590494335349E-2);
    qcopy(modelPtr_->getJointId("left_achilles_rod") - 1) = asin(t2*(-1.587289102219775E-3)+t3*1.001736520672708+t5*3.407131510426615E-4-t6*9.646678264174077E-4+1.539911055129276E-3);
    qcopy(modelPtr_->getJointId("left_ach2") - 1) = asin(t2*(-7.197863326651638E-2)+t3*8.929579539517018E-3+t5*2.669904889661342E-4-t6*8.465713056222183E-5+7.189649490089099E-2);

    // asin and acos might return nan values if input is not within a correct range
    // return error if nan value is found
    for (int i = 0; i < modelPtr_->nq; i++) {
        if (isnan(qcopy(i))) {
            throw std::runtime_error("nan values found in setupJointPositions!");
        }
    }

    fk(stance_foot_T, *model, jtype, contact_joint_id, modelPtr_->getJointId("Rz"), qcopy, stance_foot_endT, startT);
    dynamicsConstraintsPtr_->fkhofPtr_->fk(stance_foot_T, *modelPtr_, jtype, contact_joint_id, 0, q, stance_foot_endT, startT);
    Transform torso_T = stance_foot_T_des * stance_foot_T.inverse();
    qcopy.block(0, 0, 6, 1) = Transform2xyzrpy(torso_T);

    // cout << qcopy.transpose() << endl;

    // for (int i = 1; i <= modelPtr_->nq; i++) {
    //     cosq[i] = cos(qcopy(i-1));
    //     sinq[i] = sin(qcopy(i-1));
    // }
    // Eigen::VectorXd c = get_c();
    // cout << c.transpose() << endl;

    // gsl multidimensional root-finding
    const gsl_multiroot_fdfsolver_type *T;
    gsl_multiroot_fdfsolver *s;

    int status;
    size_t i, iter = 0;

    const size_t n = numDependentJoints;
    gsl_multiroot_function_fdf f = {&fillDependent_f,
                                    &fillDependent_df,
                                    &fillDependent_fdf,
                                    n, this};

    gsl_vector *x = gsl_vector_alloc(n);

    for (int i = 0; i < numDependentJoints; i++) {
        gsl_vector_set(x, i, qcopy(dynamicsConstraintsPtr_->dependentJointIds[i]));
    }

    T = gsl_multiroot_fdfsolver_hybridsj;
    s = gsl_multiroot_fdfsolver_alloc(T, n);
    gsl_multiroot_fdfsolver_set(s, &f, x);

    // // numerical check
    // gsl_matrix* J = gsl_matrix_alloc(n, n);
    // gsl_vector* f1 = gsl_vector_alloc(n);
    // gsl_vector* f2 = gsl_vector_alloc(n);
    // fillDependent_df(x, this, J);
    // for (size_t i = 0; i < n; i++) {
    //     double xx = gsl_vector_get(x, i);
    //     xx -= 1e-7;
    //     gsl_vector_set(x, i, xx);
    //     fillDependent_f(x, this, f1);
    //     xx += 2e-7;
    //     gsl_vector_set(x, i, xx);
    //     fillDependent_f(x, this, f2);
    //     xx -= 1e-7;
    //     gsl_vector_set(x, i, xx);
    //     for (size_t j = 0; j < n; j++) {
    //         double diff = gsl_matrix_get(J, j, i) - (gsl_vector_get(f2, j) - gsl_vector_get(f1, j)) / (2e-7);
    //         printf("%e ", gsl_matrix_get(J, j, i) - (gsl_vector_get(f2, j) - gsl_vector_get(f1, j)) / (2e-7));
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    // gsl_matrix_free(J);
    // gsl_vector_free(f1);
    // gsl_vector_free(f2);
    // // numerical check

    // printf("GSL ROOT FINDER START!\n");

    do {
        iter++;
        status = gsl_multiroot_fdfsolver_iterate(s);

        // double total_f = 0;
        // printf("%ld ", iter);
        // for (int i = 0; i < numDependentJoints; i++) {
        //     total_f += abs(gsl_vector_get(s->f, i));    
            // printf("%e ", abs(gsl_vector_get(s->f, i)));
        // }
        // printf("\n");
        // for (int i = 0; i < numDependentJoints; i++) {
        //     printf("%e ", abs(gsl_vector_get(s->x, i)));
        // }
        // printf("\n");
        // printf("%ld, %e\n", iter, total_f);
        // printf("\n");

        if (status) break;

        status = gsl_multiroot_test_residual(s->f, 1e-14);
    }
    while (status == GSL_CONTINUE && iter < 1000);

    // printf ("total iter = %ld, status = %s\n", iter, gsl_strerror(status));

    // at this moment, qcopy, cosq, sinq have already saved 
    // the optimal solution found by gsl!

    for (int i = 0; i < numDependentJoints; i++) {
        qcopy(dynamicsConstraintsPtr_->dependentJointIds[i]) = gsl_vector_get(s->x, i);
    }

    // gsl might be numerically unstable for some very edge cases
    // return error if found anything weird
    for (int i = 0; i < modelPtr_->nq; i++) {
        if (isnan(qcopy(i))) {
            // throw -1;
            // cout << "\tWarning: nan values found in datacomputation!\n";
        }
    }

    gsl_multiroot_fdfsolver_free(s);
    gsl_vector_free(x);

    q = qcopy;
}

int fillDependent_f(const gsl_vector* x, void *params, gsl_vector* f) {
    DigitDynamicsConstraints* constraintsData = (DigitDynamicsConstraints*)params;

    for (int i = 0; i < numDependentJoints; i++) {
        constraintsData->qcopy(constraintsData->dependentJointIds[i]) = gsl_vector_get(x, i);
    }

    // Matlab sxmbolic toolbox generates 1-index code
    // so here we shift from 0-index to 1-index
    for (int i = 1; i <= constraintsData->model->nq; i++) {
        constraintsData->cosq[i] = cos(constraintsData->qcopy(i-1));
        constraintsData->sinq[i] = sin(constraintsData->qcopy(i-1));
    }

    Eigen::VectorXd c = constraintsData->get_c();

    for (int i = 0; i < numDependentJoints; i++) {
        gsl_vector_set(f, i, c(i));
    }

    return GSL_SUCCESS;
}

int fillDependent_df(const gsl_vector* x, void *params, gsl_matrix* J) {
    DigitDynamicsConstraints* constraintsData = (DigitDynamicsConstraints*)params;

    for (int i = 0; i < numDependentJoints; i++) {
        constraintsData->qcopy(constraintsData->dependentJointIds[i]) = gsl_vector_get(x, i);
    }

    // Matlab sxmbolic toolbox generates 1-index code
    // so here we shift from 0-index to 1-index
    for (int i = 1; i <= constraintsData->model->nq; i++) {
        constraintsData->cosq[i] = cos(constraintsData->qcopy(i-1));
        constraintsData->sinq[i] = sin(constraintsData->qcopy(i-1));
    }

    Eigen::MatrixXd Jfill(numDependentJoints, constraintsData->model->nq);
    constraintsData->get_J(Jfill);

    for (int i = 0; i < numDependentJoints; i++) {
        for (int j = 0; j < numDependentJoints; j++) {
            gsl_matrix_set(J, i, j, Jfill(i, constraintsData->dependentJointIds[j]));
        }
    }

    return GSL_SUCCESS;
}

int fillDependent_fdf(const gsl_vector* x, void *params, gsl_vector* f, gsl_matrix* J) {
    fillDependent_f(x, params, f);
    fillDependent_df(x, params, J);

    return GSL_SUCCESS;
}

}; // namespace Digit
}; // namespace IDTO