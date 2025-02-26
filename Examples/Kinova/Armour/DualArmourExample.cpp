#include "DualArmourOptimizer.h"

using namespace RAPTOR;
using namespace Kinova;
using namespace Armour;

int main() {
    std::srand(std::time(nullptr));

    #ifdef NUM_THREADS
        omp_set_num_threads(NUM_THREADS);
    #else
        throw std::runtime_error("macro NUM_THREADS is not defined!");
    #endif

// FIRST ROBOT INITIALIZATION
    // read robot model and info
    const std::string robot_model_file1 = "../Robots/kinova-gen3/gen3_2f85_fixed.urdf";
    const std::string robot_info_file1 = "../Examples/Kinova/Armour/KinovaWithGripperInfo.yaml";
    const std::shared_ptr<RobotInfo> robotInfoPtr1_ = 
        std::make_shared<RobotInfo>(robot_model_file1, robot_info_file1);

    // create a trajectory instance (compute trajectory on continuous time intervals)
    Eigen::VectorXd q0 = Eigen::VectorXd::Zero(robotInfoPtr1_->num_motors);
    Eigen::VectorXd q_d0 = Eigen::VectorXd::Zero(robotInfoPtr1_->num_motors);
    Eigen::VectorXd q_dd0 = Eigen::VectorXd::Zero(robotInfoPtr1_->num_motors);

    q0 << 0.59164682,  1.05519398,  2.46569617, -1.02732907,  0.60523465,  1.11747659, -1.55976632;

        // trajectory parameters and their ranges
    Eigen::VectorXd k_center = Eigen::VectorXd::Zero(robotInfoPtr1_->num_motors);
    Eigen::VectorXd k_range = M_PI / 24 * Eigen::VectorXd::Ones(robotInfoPtr1_->num_motors);

        // trajectory duration
    const double duration = 3.0;

    std::shared_ptr<BezierCurveInterval> trajPtr1_ = 
        std::make_shared<BezierCurveInterval>(
            q0, q_d0, q_dd0, 
            k_center, k_range, 
            duration, 
            robotInfoPtr1_);
    
    // create a PZDynamics instance to compute link PZs and torque PZs
    std::shared_ptr<PZDynamics> dynPtr1_ = 
        std::make_shared<PZDynamics>(robotInfoPtr1_, trajPtr1_);

// SECOND ROBOT INITIALIZATION
    // read robot model and info
    const std::string robot_model_file2 = "../Robots/kinova-gen3/gen3_2f85_fixed_the_other_side.urdf";
    const std::string robot_info_file2 = "../Examples/Kinova/Armour/KinovaWithGripperInfo.yaml";
    const std::shared_ptr<RobotInfo> robotInfoPtr2_ = 
        std::make_shared<RobotInfo>(robot_model_file2, robot_info_file2);

    // create a trajectory instance (compute trajectory on continuous time intervals)
    q0 = Eigen::VectorXd::Ones(robotInfoPtr2_->num_motors);
    q_d0 = Eigen::VectorXd::Zero(robotInfoPtr2_->num_motors);
    q_dd0 = Eigen::VectorXd::Zero(robotInfoPtr2_->num_motors);

    q0 << 0.59164682,  1.05519398,  2.46569617, -1.02732907,  0.60523465,  1.11747659, -1.55976632;

        // trajectory parameters and their ranges
    k_center = Eigen::VectorXd::Zero(robotInfoPtr2_->num_motors);
    k_range = M_PI / 24 * Eigen::VectorXd::Ones(robotInfoPtr2_->num_motors);

    std::shared_ptr<BezierCurveInterval> trajPtr2_ =
        std::make_shared<BezierCurveInterval>(
            q0, q_d0, q_dd0, 
            k_center, k_range, 
            duration, 
            robotInfoPtr2_);

    // create a PZDynamics instance to compute link PZs and torque PZs
    std::shared_ptr<PZDynamics> dynPtr2_ = 
        std::make_shared<PZDynamics>(robotInfoPtr2_, trajPtr2_);

    // define obstacles
    const int num_obstacles = 5;
    std::vector<Eigen::Vector3d> boxCenters;
    std::vector<Eigen::Vector3d> boxOrientation;
    std::vector<Eigen::Vector3d> boxSize;

    boxCenters.resize(num_obstacles);
    boxOrientation.resize(num_obstacles);
    boxSize.resize(num_obstacles);

    for (int i = 0; i < num_obstacles; i++) {
        boxCenters[i] << 10, 0, 0.3 * i;
        boxOrientation[i] << 0, 0, 0;
        boxSize[i] << 0.1, 0.1, 0.1;
    }

    // targets
    Eigen::VectorXd q_des = Eigen::VectorXd::Random(
        robotInfoPtr1_->num_motors + robotInfoPtr2_->num_motors);

    q_des.head(robotInfoPtr1_->num_motors) << 0.05490901,  1.18798053,  2.03709998, -1.01553996, -1.45212933,  0.94367354, 0.46641401;
    q_des.tail(robotInfoPtr2_->num_motors) << 0.05490901,  1.18798053,  2.03709998, -1.01553996, -1.45212933,  0.94367354, 0.46641401;

    double t_plan = 0.5;

// COMPUTATION IN ARMOUR 
    // generate Joint Trajectory Reachable Sets
    auto start1 = std::chrono::high_resolution_clock::now();
    dynPtr1_->compute();
    dynPtr2_->compute();
    auto end1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken to generate reachable sets: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count() 
              << " ms" << std::endl;

// OPTIMIZATION
    SmartPtr<DualArmourOptimizer> mynlp = new DualArmourOptimizer();
    try {
	    mynlp->set_parameters(q_des, 
                              t_plan, 
                              robotInfoPtr1_, 
                              trajPtr1_, 
                              dynPtr1_, 
                              robotInfoPtr2_,
                              trajPtr2_,
                              dynPtr2_,
                              boxCenters,
                              boxOrientation,
                              boxSize);
    }
    catch (int errorCode) {
        throw std::runtime_error("Error initializing Ipopt! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-5);
	app->Options()->SetNumericValue("max_wall_time", 10.0);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", "ma57");
    app->Options()->SetStringValue("ma57_automatic_scaling", "yes");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    // // For gradient checking
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-8);
    // app->Options()->SetNumericValue("derivative_test_tol", 5e-4);

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded ) {
        throw std::runtime_error("Error during initialization!");
    }

    try {
        // Ask Ipopt to solve the problem
        status = app->OptimizeTNLP(mynlp);
    }
    catch (int errorCode) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    std::cout << mynlp->solution << std::endl;

    return 0;
}