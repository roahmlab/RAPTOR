#include "EndEffectorParametersIdentification.h"

using namespace RAPTOR;

const std::string folder_name = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_effector_params_data/";

Eigen::MatrixXd theta_to_LMI(Eigen::VectorXd& theta);
double BuresWassersteinDistance(Eigen::MatrixXd& A, Eigen::MatrixXd& B);
Eigen::MatrixXd matrixSquareRoot(Eigen::MatrixXd& A);
Eigen::VectorXd x_to_theta(Eigen::VectorXd& x);

int main(int argc, char* argv[]) {
    // macro NUM_THREADS should be define in cmake
    #ifdef NUM_THREADS
        omp_set_num_threads(NUM_THREADS);
    #else
        throw std::runtime_error("macro NUM_THREADS is not defined!");
    #endif

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
    const std::string urdf_filename  = "/workspaces/RAPTOR/Robots/kinova-gen3/kinova.urdf";

    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);

    Eigen::VectorXd phi = Eigen::VectorXd::Zero(10 * model.nv);
    for (int i = 0; i < model.nv; i++) {
        const int pinocchio_joint_id = i + 1;
        phi.segment<10>(10 * i) =
            model.inertias[pinocchio_joint_id]
                .toDynamicParameters();
    }

    // load the data
    std::string trajectory_filename = folder_name + "2024_11_17_no_gripper_id_" + std::string(argv[1]) + ".txt";
    // std::string trajectory_filename = "../Examples/Kinova/SystemIdentification/ExcitingTrajectories/data/T10_d5_slower/exciting-trajectory-" + std::string(argv[1]) + ".csv";

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
    sensor_noise.acceleration_error.setConstant(0.05);
  
    // Initialize the Ipopt problem
    SmartPtr<EndEffectorParametersIdentification> mynlp = new EndEffectorParametersIdentification();
    try {
	    mynlp->set_parameters(model, 
                              trajectory_filename,
                              sensor_noise,
                              H,
                              downsample_rate,
                              offset);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-8);
	app->Options()->SetNumericValue("max_wall_time", 100);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetIntegerValue("max_iter", 100);
    app->Options()->SetStringValue("mu_strategy", "monotone");
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

    std::cout <<"solution x" << mynlp->solution.transpose()<< std::endl;
    std::cout << "paramer solution: " << x_to_theta(mynlp->solution).transpose() << std::endl;
    
    Eigen::VectorXd phi_tail = phi.tail(10);
    Eigen::MatrixXd LMI_original = theta_to_LMI(phi_tail);
    Eigen::VectorXd theta_sol = x_to_theta(mynlp->solution);
    Eigen::MatrixXd LMI_solution = theta_to_LMI(theta_sol);
    double distance = BuresWassersteinDistance(LMI_original, LMI_solution);
    std::cout << "Bures-Wasserstein distance: " << distance << std::endl;

    // Utils::writeEigenMatrixToFile(mynlp->A, folder_name + "A.csv");
    // Utils::writeEigenMatrixToFile(mynlp->b, folder_name + "b.csv");

    return 0;
}

// helfer functions
double BuresWassersteinDistance(Eigen::MatrixXd& A, Eigen::MatrixXd& B) {
    if (A.rows() != A.cols() || B.rows() != B.cols() || A.rows() != B.rows()) {
        throw std::invalid_argument("Matrices A and B must be square and of the same size.");
    }

    Eigen::MatrixXd sqrtA = matrixSquareRoot(A);
    Eigen::MatrixXd middle = sqrtA * B * sqrtA;
    Eigen::MatrixXd sqrtMiddle = matrixSquareRoot(middle);

    double traceA = A.trace();
    double traceB = B.trace();
    double traceSqrtMiddle = sqrtMiddle.trace();

    return std::sqrt(traceA + traceB - 2 * traceSqrtMiddle);
}

Eigen::MatrixXd matrixSquareRoot(Eigen::MatrixXd& matrix) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();

    // Compute the square root of the diagonal matrix of eigenvalues
    Eigen::VectorXd sqrtEigenvalues = eigenvalues.array().sqrt();

    // Reconstruct the square root matrix
    return eigenvectors * sqrtEigenvalues.asDiagonal() * eigenvectors.transpose();
}

Eigen::MatrixXd theta_to_LMI(Eigen::VectorXd& theta) {
    // Extract mass
    double mass = theta(0);

    // Extract center of mass (com)
    Eigen::Vector3d com = theta.segment<3>(1);

    // Extract inertia elements
    double Ixx = theta(4);
    double Ixy = theta(5);
    double Iyy = theta(6);
    double Ixz = theta(7);
    double Iyz = theta(8);
    double Izz = theta(9);

    // Construct the inertia matrix
    Eigen::Matrix3d inertia;
    inertia << Ixx, Ixy, Ixz,
               Ixy, Iyy, Iyz,
               Ixz, Iyz, Izz;

    // Compute trace of the inertia matrix
    double inertia_trace = inertia.trace();

    // Compute the LMI matrix
    Eigen::MatrixXd LMI(4, 4);
    LMI.setZero();
    LMI.block<3, 3>(0, 0) = 0.5 * inertia_trace * Eigen::Matrix3d::Identity() - inertia;
    LMI.block<3, 1>(0, 3) = com;                                                        
    LMI.block<1, 3>(3, 0) = com.transpose();                                          
    LMI(3, 3) = mass;                                                                  

    return LMI;
}

Eigen::VectorXd x_to_theta(Eigen::VectorXd& z) {
    double d1 = z[0];
    double d2 = z[1];
    double d3 = z[2];
    double d4 = z[3];
    double s12 = z[4];
    double s23 = z[5];
    double s13 = z[6];
    double t1 = z[7];
    double t2 = z[8];
    double t3 = z[9];

    Eigen::Matrix4d U;
    U << std::exp(d1), s12,      s13,      t1,
        0.0,     std::exp(d2),  s23,      t2,
        0.0,     0.0,      std::exp(d3),  t3,
        0.0,     0.0,      0.0,      std::exp(d4);

    // Compute LMI = U' * U
    Eigen::Matrix4d LMI = U.transpose() * U;

    // End-effector parameters
    Eigen::VectorXd theta = Eigen::VectorXd::Zero(10);
    theta(0) = LMI(3, 3);                        
    theta.segment<3>(1) = LMI.block<3, 1>(0, 3); 
    theta(4) = LMI(1, 1) + LMI(2, 2);            // IXX
    theta(5) = -LMI(0, 1);                       // IXY
    theta(6) = LMI(0, 0) + LMI(2, 2);            // IYY
    theta(7) = -LMI(0, 2);                       // IXZ
    theta(8) = -LMI(1, 2);                       // IYZ
    theta(9) = LMI(0, 0) + LMI(1, 1);            // IZZ

    return theta;
}
