#include "SafePayloadExcitingTrajectoryGenerator.h"

using namespace RAPTOR;
using namespace Kinova;
using namespace Armour;
using namespace Ipopt;

int main(int argc, char* argv[]) {
    #ifdef NUM_THREADS
        omp_set_num_threads(NUM_THREADS);
    #else
        throw std::runtime_error("macro NUM_THREADS is not defined!");
    #endif

    const bool use_momentum_regressor_or_not = true;

    // INITIALIZATION
    // read robot model and info
    const std::string robot_model_file = "../Robots/kinova-gen3/kinova_grasp_fixed.urdf";
    const std::string robot_info_file = "../Examples/Kinova/Armour/KinovaWithGripperInfo.yaml";
    const std::shared_ptr<RobotInfo> robotInfoPtr_ = 
        std::make_shared<RobotInfo>(robot_model_file, robot_info_file);

    // create a trajectory instance (compute trajectory on continuous time intervals)
        // initial conditions of the trajectory
    const Eigen::VectorXd q0 = Eigen::VectorXd::Random(robotInfoPtr_->num_motors);
    const Eigen::VectorXd q_d0 = Eigen::VectorXd::Random(robotInfoPtr_->num_motors);
    const Eigen::VectorXd q_dd0 = Eigen::VectorXd::Random(robotInfoPtr_->num_motors);

        // trajectory parameters and their ranges
    const Eigen::VectorXd k_center = Eigen::VectorXd::Random(robotInfoPtr_->num_motors);
    const Eigen::VectorXd k_range = M_PI / 24 * Eigen::VectorXd::Ones(robotInfoPtr_->num_motors);

        // trajectory duration
    const double duration = 3.0;

    std::shared_ptr<BezierCurveInterval> trajIntervalPtr_ = 
        std::make_shared<BezierCurveInterval>(
            q0, q_d0, q_dd0, 
            k_center, k_range, 
            duration, 
            robotInfoPtr_);
    
    // create a PZDynamics instance to compute link PZs and torque PZs
    std::shared_ptr<PZDynamics> dynPtr_ = 
        std::make_shared<PZDynamics>(robotInfoPtr_, trajIntervalPtr_);

    std::vector<Eigen::Vector3d> boxCenters = {
        Eigen::Vector3d(0.0, 0.0, 0.15), // ground
        Eigen::Vector3d(0.53, 0.49, 0.56),  // back wall
        Eigen::Vector3d(-0.39, -0.84, 0.56), // bar near the control
        Eigen::Vector3d(-0.39, -0.17, 0.56), //bar bewteen 10 and 20 change to wall
        Eigen::Vector3d(0.0, 0.0, 1.12), //ceiling
        Eigen::Vector3d(0.47, -0.09, 1.04) // top camera  
    };    
    std::vector<Eigen::Vector3d> boxOrientations = {
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0)
    };
    std::vector<Eigen::Vector3d> boxSizes = {
        Eigen::Vector3d(5.0, 5.0, 0.01),
        Eigen::Vector3d(5.0, 0.05, 1.12),
        Eigen::Vector3d( 0.05, 0.05, 1.12),
        Eigen::Vector3d( 0.05, 1.28, 1.28),
        Eigen::Vector3d( 5, 5, 0.05),
        Eigen::Vector3d( 0.15, 0.15, 0.15)
    };   

// COMPUTATION IN ARMOUR 
    // generate Joint Trajectory Reachable Sets
    auto start1 = std::chrono::high_resolution_clock::now();
    dynPtr_->compute();
    auto end1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken to generate reachable sets: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count() 
              << " ms" << std::endl;

    // Initialize optimizer
    SmartPtr<SafePayloadExcitingTrajectoryGenerator> mynlp = new SafePayloadExcitingTrajectoryGenerator();
    try {
	    mynlp->set_parameters(robotInfoPtr_,
                              trajIntervalPtr_,
                              dynPtr_,
                              boxCenters,
                              boxOrientations,
                              boxSizes,
                              use_momentum_regressor_or_not);
        mynlp->constr_viol_tol = 1e-5;
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-6);
    app->Options()->SetNumericValue("constr_viol_tol", mynlp->constr_viol_tol);
	app->Options()->SetNumericValue("max_wall_time", 5.0);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", "ma57");
    app->Options()->SetStringValue("ma57_automatic_scaling", "yes");
	if (mynlp->enable_hessian) {
        app->Options()->SetStringValue("hessian_approximation", "exact");
    }
    else {
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    }

    // // For gradient checking
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-3);

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded ) {
		throw std::runtime_error("Error during initialization of optimization!");
    }

    // Run ipopt to solve the optimization problem
    double solve_time = 0;
    try {
        auto start = std::chrono::high_resolution_clock::now();

        // Ask Ipopt to solve the problem
        status = app->OptimizeTNLP(mynlp);

        auto end = std::chrono::high_resolution_clock::now();
        solve_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Total solve time: " << solve_time << " milliseconds.\n";
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    return 0;
}