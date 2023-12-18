#include "DigitSingleStepOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include <iomanip>

using namespace IDTO;
using namespace Digit;
using namespace Ipopt;

using std::cout;
using std::endl;

int main() {
    // Eigen::initParallel();

    // define robot model
    const std::string urdf_filename = "../Examples/Digit/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = -9.806;

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    Eigen::VectorXi jtype(model.nq);
    jtype << 4, 5, 6, 1, 2, 3, 
             3, 3, -3, 3, 2, 3, 3, 3, 3, 2, 3, 3, 2, 3, 3,
             3, 3, -3, 3, 2, 3, 3, 3, 3, 2, 3, 3, 2, 3, 3;
    
    // ignore friction for now
    model.friction.setZero();

    // manually import motor inertia 
    model.rotorInertia(model.getJointId("left_hip_roll") - 1) = 0.173823936;
    model.rotorInertia(model.getJointId("left_hip_yaw") - 1) = 0.067899975;
    model.rotorInertia(model.getJointId("left_hip_pitch") - 1) = 0.1204731904;
    model.rotorInertia(model.getJointId("left_knee") - 1) = 0.1204731904;
    model.rotorInertia(model.getJointId("left_toe_A") - 1) = 0.036089475;
    model.rotorInertia(model.getJointId("left_toe_B") - 1) = 0.036089475;
    model.rotorInertia(model.getJointId("right_hip_roll") - 1) = 0.173823936;
    model.rotorInertia(model.getJointId("right_hip_yaw") - 1) = 0.067899975;
    model.rotorInertia(model.getJointId("right_hip_pitch") - 1) = 0.1204731904;
    model.rotorInertia(model.getJointId("right_knee") - 1) = 0.1204731904;
    model.rotorInertia(model.getJointId("right_toe_A") - 1) = 0.036089475;
    model.rotorInertia(model.getJointId("right_toe_B") - 1) = 0.036089475;

    pinocchio::Data data(model);

    const double T = 0.4;
    const int N = 48;
    const int degree = 6;

    GaitParameters gp;
    gp.swingfoot_midstep_z_des = 0.30;
    gp.swingfoot_begin_y_des = 0.40;
    gp.swingfoot_end_y_des = -0.40;

    FourierCurves fc(T, N, NUM_INDEPENDENT_JOINTS, Chebyshev, degree);

    std::ifstream initial_guess("initial-digit.txt");

    Eigen::VectorXd z(252);
    for (int i = 0; i < 252; i++) {
        initial_guess >> z(i);
    }
    initial_guess.close();
    // z << 0.29224492320581330063, 0.84522753120149440864, -0.53472035519430305417, 2.67516301546460333327, -0.44634944057612102597, -0.76557708099529064860, -0.87947584624992791724, 1.18856119611949040404, 0.67883259766877379615, 1.03747494847933618978, -1.54471213745385416161, -0.28709462548397673975, 
    //     -0.22619963156456529552, 8.32970255735735598535, 0.11305008784314241055, 0.60920762849420551088, -0.69267911761540668092, 1.56258951342091512160, 0.56788975731392854485, 0.04509102232474253985, -0.51504998170545313041, -0.26225866224635269175, 0.60807085247185110877, -0.31384258952058263770, 
    //     -0.50575334876664701511, 0.65948586832083855302, -0.69864503520547860393, 9.60939323797039257613, 0.39340895915250284620, 2.11979974666413184536, -2.94739786631211853063, -0.55759087450089106497, 1.63810229809724683392, 0.27065783608553611872, 2.30217582696675737708, 0.54400235272765629091, 
    //     2.57075404598062240424, 0.30464193289494667649, 1.09194193028561969960, -0.01740614979035023241, 0.07566677625570045196, 12.11110988046748992986, -1.39050180540355672321, -1.98808371012529372379, 2.29262962401298597825, -1.02867039839934060197, 0.68175614088205160890, 0.24024103707123478646, 
    //     0.35837565969287449308, -0.97372325285554639507, 1.73706602245833763298, 0.69314722550882545971, 1.94701024167037117785, 0.65653030083032715147, 0.58907499641362859322, 10.86991224084809637418, -0.00972739272540661636, 0.54326957673493303513, 0.15719870761440069007, 0.95928540245879545534, 
    //     -0.17719247715112562691, 0.35116358324194901774, -0.71977736084598276101, 0.23351365839646306966, -1.46514843349017787055, -0.44060915941917350169, -1.28423685808827325872, -0.32007628208818605930, -0.35912343773621857590, 10.81049519027172323149, -0.93966340177294727098, -1.12396398626085392003, 
    //     0.75533724389428658075, -0.64325738087021844258, 0.95651408640387136639, 0.41433726443426260877, 0.80420847517750293498, 0.12435011911263518791, 1.06813677211686552937, 0.34670507142461520766, 1.07170818714559956675, 0.25081148319876644948, 0.34344421588577395132, 10.85807062524941279946, 
    //     0.11109146887103597823, -0.84834384226477022040, 0.66541026041380646472, -0.58263811498498840891, -2.48714623151437175252, -0.19220565717908852377, -0.34687497932001681855, -0.16172991979402998042, 0.37635683762946570141, 0.11256661300020291694, -0.02717458499107647252, 0.02100945014581665979, 
    //     0.01742831194261194871, 11.76901098259756928144, -0.31156918315494569471, -1.03467727606208192981, 0.62328975878576509118, -0.07227534563290326231, -0.01012906449374933321, 0.00426236564982853446, 0.55668369923494298579, -0.80448196056949949906, -1.24554887494011023996, 0.41017830569240676386, 
    //     0.04629673726445340115, 0.01878592757874619421, -0.01738054810634838757, 12.57301274264597878982, 8.97351257773206967272, 26.87040885282691249358, -27.66190233475770554605, 8.19233197125679168948, 12.63508088910562143781, 4.68065175902507846217, 0.53105208067470599520, -1.44939784103564850426, 
    //     -2.02567814693920933422, -0.62936782397591195526, -0.44795755376786372537, 0.03725091644627713527, 0.12357215310273643449, 11.59003495815929518642, 11.25772250628561721442, 51.71653544170943206382, -45.35022117126981555657, 6.28137256869538074255, 26.40974813850033342533, 7.17325462266716407100, 
    //     11.14028337699130943861, 2.99378405510431910841, 1.44017901450531748964, 1.32203180022382027481, 0.69810338128694648940, 0.00441416615956360427, 1.04747526105167532151, 12.25593088605299918470, -5.29884244178106822432, -2.80297406982991503810, 2.74145171708251211840, -25.61268126428158709018, 
    //     18.83866881892104672147, -17.86032627970400454842, -6.02206322427582918522, -0.20664561366472117276, -1.40480858791497298377, -3.12025270039467450545, -0.84749487265284362270, -1.26369220799312320658, -2.09113209623398610404, 6.32698577531110917249, 5.32386034373057892566, 4.08032178976901871437, 
    //     -6.66181968818669023591, 29.86401910237904999690, -15.43452753185406400860, 7.38704331467935038091, 7.52363576985012194598, 1.82033381545548844649, 0.10072027455891263692, 2.20964651762191177653, 1.97123704022434997896, -0.26843866986070147318, 1.40225471018570790882, 7.04644582771768579477, 
    //     0.35656014637352806407, -0.00906395854354049603, 0.46069946696759167537, 0.74839096529395809920, -0.33533878513553261325, 0.32128648704649148682, -0.35823383649525158967, 0.01277719633451589870, -0.46143196171624611956, -0.74898843262303371748, 0.33597173590593826242, -0.32203949780971108385, 
    //     0.04048147829951193738, 0.01536230567254849577, -0.01140340707413285345, -0.03624201348997058841, 0.05196481209930344652, 0.02772111307592835250, 0.04123346241206938406, 0.01575567147071101765, -0.01332605035037683419, -0.03710505620828823126, 0.05272215286200387929, 0.02751321417440268438, 
    //     -0.10757221480150275239, -0.00008915734807600152, 0.00039994304199106361, -0.00050045587234690523, -0.05974272299413958126, 0.04597751540479672050, -0.04123710838115582367, -0.01575568183317761425, 0.01332759339521411360, 0.03721645734077697054, 0.00204089428979350399, 0.03711188024704405869, 
    //     -0.03790098120502969486, -0.05272204248532023213, 0.05537759865569773232, 0.00164559913531180601, -0.02751319005977243340, 0.02944450274572157217, -0.00529187612345230051, 0.01301355086131937393, -0.11823255591477874193, -0.04048036941580287346, -0.01536015467525131926, 0.01140283804432499354, 
    //     0.03634414176274806257, -0.00199214717740287027, 0.03624010836234700661, -0.03702100599641015150, -0.05195571906107050025, 0.05458271341568621260, -0.00165066837248612906, -0.02771228007229493973, 0.02964061353463749368, 0.00527764831518921958, 0.01252063369794281615, -0.11750211994409595628, 
    //     -0.02475375104448079497, -0.00076026424322239087, 0.00085584228616659161, -0.01183458399773307593, -0.00052183382242302679, 0.00092535234400804162, 0.00050008464763759452, 0.00001564471435192455, -0.00002485827474769738, -0.00046179753123938577, 0.00002254106914544711, -0.00002361026544127246, 
    //     0.07756550781211253576, 0.06403436178635335541, -0.00652352390272287978, -0.00112674019614598334, 0.00090440306802109550, 0.00016222614670252254, 0.05192248609980511243, -0.00417357303760403000, 0.06306215728958201872, 0.00210241566465715207, 0.00157469136643473458, -0.00262199004698662926;

    SmartPtr<DigitSingleStepOptimizer> mynlp = new DigitSingleStepOptimizer();
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
	app->Options()->SetNumericValue("max_wall_time", 500);
    app->Options()->SetNumericValue("obj_scaling_factor", 1e-5);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-4);
    app->Options()->SetIntegerValue("max_iter", 200);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetStringValue("mu_strategy", "monotone");
    app->Options()->SetStringValue("linear_solver", "ma57");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetStringValue("nlp_scaling_method", "none");

    // For gradient checking
    app->Options()->SetStringValue("output_file", "ipopt_digit.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("point_perturbation_radius", 1e-3);
    // // app->Options()->SetIntegerValue("derivative_test_first_index", 168);
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-4);

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded ) {
		throw std::runtime_error("Error during initialization of optimization!");
    }

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
    
    // Index n, m, nnz_jac_g, nnz_h_lag;
    // TNLP::IndexStyleEnum index_style;
    // mynlp->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
    // Number g[mynlp->numCons];
    // Number values[mynlp->numCons * mynlp->numVars];
    // // mynlp->eval_g(mynlp->numVars, x, false, mynlp->numCons, g);
    // mynlp->eval_jac_g(mynlp->numVars, x, false, mynlp->numCons, 0, NULL, NULL, values);

    // Print the solution
    if (mynlp->solution.size() == mynlp->numVars) {
        std::ofstream solution("solution-digit.txt");
        solution << std::setprecision(20);
        for (int i = 0; i < mynlp->numVars; i++) {
            solution << mynlp->solution[i] << std::endl;
        }
        solution.close();

        std::ofstream trajectory("trajectory-digit.txt");
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

    return 0;
}