#include "EndEffectorParametersIdentificationMomentum.h"

using namespace RAPTOR;

const std::string folder_name = "../Examples/Kinova/SystemIdentification/ParametersIdentification/end_effector_params_data/";

int main(int argc, char* argv[]) {
    // check if the file number is provided
    if (argc < 2) {
        throw std::invalid_argument("No arguments provided. please choose the trajectory file");
    }

    std::string file_number = std::string(argv[1]);
    int H = 5;

    if (argc > 2) {
        H = std::stoi(argv[2]);
    }

    int downsample_rate = 1;

    if (argc > 3) {
        downsample_rate = std::stoi(argv[3]);
    }

    // Load the robot model
    const std::string urdf_filename  = "../Robots/kinova-gen3/gen3_2f85_fixed_object.urdf";

    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);
    
    // load the data
    std::string trajectory_filename = folder_name + "2024_11_17_no_gripper_id_" + std::string(argv[1]) + ".txt";
    // std::string trajectory_filename = "../Examples/Kinova/SystemIdentification/ExcitingTrajectories/data/T10_d5_slower/exciting-trajectory-" + std::string(argv[1]) + ".csv";
    // std::string trajectory_filename = "../torque_output.txt";

    // load friction parameters
    const std::string friction_parameters_filename = folder_name + "friction_params.csv";
    Eigen::VectorXd friction_parameters = Utils::initializeEigenMatrixFromFile(friction_parameters_filename).col(0);
    if (friction_parameters.size() != 3 * model.nv &&
        friction_parameters.size() != 4 * model.nv) {
        std::cerr << "Friction_parameters size: " << friction_parameters.size() << std::endl;
        throw std::runtime_error("Friction solution file is wrong");
    }

    model.friction = friction_parameters.head(model.nv);
    model.damping = friction_parameters.segment(model.nv, model.nv);
    model.armature  = friction_parameters.segment(2 * model.nv, model.nv);
    Eigen::VectorXd offset = Eigen::VectorXd::Zero(model.nv);
    if (friction_parameters.size() == 4 * model.nv) {
        offset = friction_parameters.tail(model.nv);
    }

    // Sensor noise info
    SensorNoiseInfo sensor_noise(model.nv);
    // sensor_noise.position_error << Utils::deg2rad(0.02),
    //                                Utils::deg2rad(0.02),
    //                                Utils::deg2rad(0.02),
    //                                Utils::deg2rad(0.02),
    //                                Utils::deg2rad(0.011),
    //                                Utils::deg2rad(0.011),
    //                                Utils::deg2rad(0.011); // from Kinova official support, resolution of joint encoders
    // sensor_noise.velocity_error = 5 * sensor_noise.position_error;
    sensor_noise.position_error.setZero();
    sensor_noise.velocity_error.setZero();
    sensor_noise.acceleration_error_type = SensorNoiseInfo::SensorNoiseType::Ratio;
    sensor_noise.acceleration_error.setConstant(0.10);
  
    // Initialize the Ipopt problem
    SmartPtr<EndEffectorParametersIdentificationMomentum> mynlp = new EndEffectorParametersIdentificationMomentum();

    double setup_time = 0;
    try {
        auto start = std::chrono::high_resolution_clock::now();
	    mynlp->set_parameters(model,
                              offset);
        mynlp->add_trajectory_file(trajectory_filename,
                                   sensor_noise,
                                   H,
                                   TimeFormat::Nanosecond,
                                   downsample_rate);
        auto end = std::chrono::high_resolution_clock::now();
        setup_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Setup time: " << setup_time << " milliseconds.\n";
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-10);
	app->Options()->SetNumericValue("max_wall_time", 10.0);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetIntegerValue("max_iter", 100);
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
    // app->Options()->SetStringValue("derivative_test", "second-order");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-8);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-5);
    // app->Options()->SetNumericValue("point_perturbation_radius", 1);

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
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    std::cout << "solution: " << mynlp->solution.transpose()<< std::endl;
    std::cout << "parameter solution: " << mynlp->z_to_theta(mynlp->solution).transpose() << std::endl;
    std::cout << "groundtruth: " << mynlp->phi_original.tail(10).transpose() << std::endl;

    return 0;
}
