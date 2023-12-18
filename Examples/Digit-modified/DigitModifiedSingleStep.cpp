#include "DigitModifiedSingleStepOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include <iomanip>

using namespace IDTO;
using namespace DigitModified;
using namespace Ipopt;

using std::cout;
using std::endl;

int main() {
    // Eigen::initParallel();

    // Define robot model
    const std::string urdf_filename = "../Examples/Digit-modified/digit-v3-modified.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = -9.806;

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    Eigen::VectorXi jtype(model.nq);
    jtype << 4, 5, 6, 1, 2, 3, 
             3, 3, -3, 3, 3, 3, 3,
             3, 3, -3, 3, 3, 3, 3;
    
    // ignore all motor dynamics
    model.rotorInertia.setZero();
    model.damping.setZero();
    model.friction.setZero();

    pinocchio::Data data(model);

    const double T = 0.4;
    const int N = 32;
    const int degree = 5;

    GaitParameters gp;
    gp.swingfoot_midstep_z_des = 0.22;
    gp.swingfoot_begin_x_des = -0.25;
    gp.swingfoot_begin_y_des = 0.40;
    gp.swingfoot_end_x_des = -0.25;
    gp.swingfoot_end_y_des = -0.40;                                  

    std::ifstream initial_guess("initial-digit-modified.txt");

    Eigen::VectorXd z(222);
    for (int i = 0; i < 222; i++) {
        initial_guess >> z(i);
    }
    initial_guess.close();
    // z << 0.44001300673447579781, 0.90150841199927378433, -1.36896462498458881818, 1.17514427366831353261, 1.49127770819168103955, 2.12381489742914641994, -0.82759129050768542868, 0.53137207451734613795, -1.67995230143878249152, -0.18866375650903835504, -0.15893353665649906370, 10.17378464326555764785, 
    //     0.24906378256190339626, 1.13372571890367868086, -1.17876565189164517200, 0.32230045632300652336, 0.91608278315458502306, 0.67928722900229043802, 0.31490334064012565074, -0.63447245508795857560, -0.54140814640940027047, 0.65895200145184085194, -0.75282951196229741520, 11.50359973237056721018, 
    //     0.22939257155599229732, 1.85268430597661759585, -2.74652666636462416960, -0.59901875070369114962, 1.94036815765457304472, 0.41091855900525037226, 2.50613343617640138561, 0.73729295561515439950, 2.61055552506960308889, 0.40808039640519155888, 1.03271893957363380956, 12.08773903143128869431, 
    //     -1.33178100930289367021, -2.95648084502303376553, 2.02834495242657864011, -1.35377054087555315220, 1.06647358911332434950, 0.46617554589011150901, 2.77487968187975786449, 1.05802250689076249301, 1.36100895521274734001, 0.35631095909228799323, 0.62219453695030346729, 12.90431915773569393480, 
    //     1.36101683384519644981, 3.01993492025489507213, -2.07233937640846077599, 1.38269061354209976678, -1.09197048504790616974, -0.47712860253501543362, -2.83387032147723738618, -1.08161947760753651693, -1.38927825482159672532, -0.36433122137585216560, -0.63491662392911463186, 12.90340824184330692503, 
    //     -1.07588550410097760768, -1.65054142078163579299, 1.18314358240267836209, -0.71606532495998920229, 1.55069614201935679354, 0.83169471446936715786, 1.83763559289506672734, 0.89741527961541456282, 1.00312768448165501312, 0.27168072541095039663, 0.38602401522200618933, 12.24724637404398208673, 
    //     -37.37400417114449169276, 11.83966118846237414175, 63.18301464776368447929, 39.26540777359542744307, -15.42102922821363897299, -10.84341721551445125726, -16.66699976802800620135, -3.96734357111079694036, 4.15012756351070155603, 0.70038338253436527303, 0.21004702104128569751, 7.14292379236334706860, 
    //     0.11037571751584196456, -0.73162763548524034629, 0.60261196581836096176, -1.06840960707211007907, -2.32563640863249343482, -0.09055876285268542325, -0.13180822091658259132, 0.16277265113199645974, 0.36535748849497046153, 0.07440960790695111404, -0.16963511610697115661, 11.28985908555694805955, 
    //     -0.39011314159654841838, -0.95479782978006366712, 0.81092153786823450989, -0.01753990148161576473, -0.02737876233333116915, 0.19303635441327296474, 0.25096116161307652215, -1.24051997525845747461, -0.62850150672428428678, 0.62218936459479667889, -0.35479363473799746531, 11.86174797500852839960, 
    //     8.48254989914047463628, 26.63144097910578267374, -26.97271222226235920516, 7.89184302863042486109, 13.26510099634745998287, 5.08270107407477400585, 1.06033670864361018182, -0.95318090798395249674, -2.12001725165167087539, -0.48815885865954405931, -0.66692167811402780497, 11.63960439553620851427, 
    //     8.83999453657376044191, 66.64794017120411240285, -36.11718148461645938596, 0.79510611592550151094, 15.91434894117091403132, 0.62018309221661105379, 10.11443860068377453842, 2.36685895569061210963, 3.25188253705358754075, 1.62099396814928975452, 1.70073192248763715995, 13.47698843902188059474, 
    //     -9.09566829502284690534, -67.26737599326115457643, 36.78424897033493579102, -1.19509592455893942642, -15.70236938060629050540, -0.82238729717296843447, -10.32273548564522691606, -2.52512466883874120427, -3.40322401621424752705, -1.72462411303768492488, -1.74926036072927693077, 13.45428343967937756531, 
    //     10.31693182960247412439, 8.53430781136323091118, -24.30963629384957869206, 31.15774929019300643063, 5.65616207138983551062, 1.28963672020117936512, 0.72190483156029383593, 1.99647528882758318325, 2.19673134655714186181, -0.57198991080049399471, 0.69827195277475861968, 8.87196979308638056239, 
    //     -2.16256563250479416638, -2.93474138002620188459, 4.25425448987740040252, -2.15180217810646823295, -2.10460587814073418755, -0.64850608024385136829, 0.36870250624374023563, -0.04276990222657374441, 0.95696739023226218279, -0.19009528847653586747, -0.10605211312320921013, 10.40477471911005658001, 
    //     0.35656014637352806407, -0.00906395854354049603, 0.46069946696759167537, 0.74839096529395809920, -0.76457344670744886717, 0.33217315300430039215, -0.02101876539890849133, -0.35823383649525158967, 0.01277719633451589870, -0.46143196171624611956, -0.74898843262303371748, 0.76518348711330030909, 
    //     -0.33285746780955938195, 0.02083713516481719083, 0.04048147829951193044, 0.01536230567254849577, -0.01140340707413284999, -0.03624201348997058147, 0.03700391780468362951, -0.01253796833328908560, 0.11746290362990911837, 0.04123346241206938406, 0.01575567147071101765, -0.01332605035037692787, 
    //     -0.03710505620828818962, 0.03788668354790423326, -0.01301446024173266193, 0.11824098014743654461, -0.10757221480150275239, -0.00008915734807600152, 0.00039994304199106361, -0.00050045587234690523, -0.05974272299413958126, 0.04597751540479672050, -0.04123710838115582367, -0.01575568183317761425, 
    //     0.01332759339521411360, 0.03711188024704405869, -0.03790098120502969486, 0.01301355086131937393, -0.11823255591477874193, -0.04048036941580287346, -0.01536015467525131926, 0.01140283804432499354, 0.03624010836234700661, -0.03702100599641015150, 0.01252063369794281615, -0.11750211994409595628, 
    //     -0.02475375104448079497, -0.00076026424322239087, 0.00085584228616659161, -0.01183458399773307593, -0.00052183382242302679, 0.00092535234400804162;

    SmartPtr<DigitModifiedSingleStepOptimizer> mynlp = new DigitModifiedSingleStepOptimizer();
    try {
	    mynlp->set_parameters(z,
                              T,
                              N,
                              degree,
                              model,
                              jtype,
                              gp);
    }
    catch (int errorCode) {
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-4);
	app->Options()->SetNumericValue("max_wall_time", 300);
    app->Options()->SetNumericValue("obj_scaling_factor", 1e-4);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-4);
    app->Options()->SetIntegerValue("max_iter", 1000);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetStringValue("mu_strategy", "monotone");
    app->Options()->SetStringValue("linear_solver", "ma86");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    app->Options()->SetStringValue("nlp_scaling_method", "none");

    // For gradient checking
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("point_perturbation_radius", 1e-2);
    // // app->Options()->SetIntegerValue("derivative_test_first_index", 168);
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-6);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-4);

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded ) {
		throw std::runtime_error("Error during initialization of optimization!");
    }

    // Run ipopt to solve the optimization problem
    try {
        auto start = std::chrono::high_resolution_clock::now();

        // Ask Ipopt to solve the problem
        status = app->OptimizeTNLP(mynlp);

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Total solve time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds.\n";
    }
    catch (int errorCode) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }
    
    // Print the solution
    if (mynlp->solution.size() == mynlp->numVars) {
        std::ofstream solution("solution-digit-modified.txt");
        solution << std::setprecision(20);
        for (int i = 0; i < mynlp->numVars; i++) {
            solution << mynlp->solution[i] << std::endl;
        }
        solution.close();

        std::ofstream trajectory("trajectory-digit-modified.txt");
        trajectory << std::setprecision(20);
        for (int i = 0; i < NUM_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->dcidPtr_->q(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->dcidPtr_->v(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->dcidPtr_->a(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->dcidPtr_->tau(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
            for (int j = 0; j < N; j++) {
                trajectory << mynlp->dcidPtr_->lambda(j)(i) << ' ';
            }
            trajectory << std::endl;
        }
        trajectory.close();
    }

    // Number x[z.size()];
    // for (int i = 0; i < z.size(); i++) {
    //     x[i] = z(i);
    // }
    // Index n, m, nnz_jac_g, nnz_h_lag;
    // TNLP::IndexStyleEnum index_style;
    // mynlp->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
    // Number f;
    // Number grad_f[mynlp->numVars];
    // mynlp->eval_f(mynlp->numVars, x, false, f);
    // mynlp->eval_grad_f(mynlp->numVars, x, false, grad_f);
    // Number g[mynlp->numCons];
    // Number values[mynlp->numCons * mynlp->numVars];
    // mynlp->eval_g(mynlp->numVars, x, false, mynlp->numCons, g);
    // mynlp->eval_jac_g(mynlp->numVars, x, false, mynlp->numCons, 0, NULL, NULL, values);

    return 0;
}