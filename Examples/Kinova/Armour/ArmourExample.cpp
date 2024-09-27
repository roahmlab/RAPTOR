#include "ArmourOptimizer.h"

using namespace RAPTOR;
using namespace Armour;

int main() {
    std::srand(std::time(nullptr));

    #ifdef NUM_THREADS
        omp_set_num_threads(NUM_THREADS);
    #else
        throw std::runtime_error("macro NUM_THREADS is not defined!");
    #endif

// INITIALIZATION
    // read robot model and info
    const std::string robot_model_file = "../Robots/kinova-gen3/kinova.urdf";
    const std::string robot_info_file = "../Examples/Kinova/Armour/KinovaWithoutGripperInfo.yaml";
    const std::shared_ptr<RobotInfo> robotInfoPtr_ = 
        std::make_shared<RobotInfo>(robot_model_file, robot_info_file);

    // create a trajectory instance (compute trajectory on continuous time intervals)
        // initial conditions of the trajectory
    const Eigen::VectorXd q0 = Eigen::VectorXd::Random(robotInfoPtr_->num_motors);
    const Eigen::VectorXd q_d0 = Eigen::VectorXd::Random(robotInfoPtr_->num_motors);
    const Eigen::VectorXd q_dd0 = Eigen::VectorXd::Random(robotInfoPtr_->num_motors);

        // trajectory parameters and their ranges
    const Eigen::VectorXd k_center = Eigen::VectorXd::Random(robotInfoPtr_->num_motors);
    const Eigen::VectorXd k_range = M_PI / 48 * Eigen::VectorXd::Ones(robotInfoPtr_->num_motors);

        // trajectory duration
    const double duration = 3.0;

    std::shared_ptr<BezierCurveInterval> trajPtr_ = 
        std::make_shared<BezierCurveInterval>(
            q0, q_d0, q_dd0, 
            k_center, k_range, 
            duration, 
            robotInfoPtr_->ultimate_bound_info);
    
    // create a KinematicsDynamics instance to compute link PZs and torque PZs
    std::shared_ptr<KinematicsDynamics> kdPtr = 
        std::make_shared<KinematicsDynamics>(robotInfoPtr_, trajPtr_);

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

    // suction cup info
    const double suction_force = 0.0; // N
    const double mu = 0.7; // friction coefficient
    const double suction_radius = 0.1; // m

    // targets
    Eigen::VectorXd q_des = Eigen::VectorXd::Random(robotInfoPtr_->num_motors);
    double t_plan = 0.5 * duration;

// COMPUTATION IN ARMOUR 
    // generate Joint Trajectory Reachable Sets
    auto start1 = std::chrono::high_resolution_clock::now();
    GenerateJRS(robotInfoPtr_, trajPtr_);
    auto end1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken to generate JRS: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count() 
              << " ms" << std::endl;
    
    // generate Link and Torque PZs
    auto start2 = std::chrono::high_resolution_clock::now();
    GenerateLinkAndTorquePZs(robotInfoPtr_, trajPtr_, kdPtr);
    auto end2 = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken to generate Link and Torque PZs: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count() 
              << " ms" << std::endl;

    // compute Robust Input Bounds
    auto start3 = std::chrono::high_resolution_clock::now();
    Eigen::MatrixXd torque_radius = ComputeRobustInputBounds(robotInfoPtr_, trajPtr_, kdPtr);
    auto end3 = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken to compute Robust Input Bounds: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count() 
              << " ms" << std::endl;

// OPTIMIZATION
    SmartPtr<ArmourOptimizer> mynlp = new ArmourOptimizer();
    try {
	    mynlp->set_parameters(q_des, 
                              t_plan, 
                              robotInfoPtr_, 
                              trajPtr_, 
                              kdPtr, 
                              torque_radius,
                              boxCenters,
                              boxOrientation,
                              boxSize,
                              suction_force,
                              mu,
                              suction_radius);
    }
    catch (int errorCode) {
        throw std::runtime_error("Error initializing Ipopt! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-5);
	app->Options()->SetNumericValue("max_wall_time", 1.0);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", "ma57");
    app->Options()->SetStringValue("ma57_automatic_scaling", "yes");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    // // For gradient checking
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-5);

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