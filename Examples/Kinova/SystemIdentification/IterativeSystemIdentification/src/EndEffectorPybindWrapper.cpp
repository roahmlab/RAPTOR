#include "EndEffectorPybindWrapper.h"

namespace RAPTOR::Kinova {

EndEffectorPybindWrapper::EndEffectorPybindWrapper(const std::string urdf_filename,
                                                   const std::string trajectory_filename,
                                                   const std::string friction_parameters_filename,
                                                   const int H_input,
                                                   const bool display_info):
                                                   H(H_input),trajectory_filename(trajectory_filename) {
    // Define robot model
    pinocchio::Model model_double;
    pinocchio::urdf::buildModel(urdf_filename, model_double);
    model = model_double.cast<double>();

    // load friction parameters
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
        std::cout << offset.transpose() << std::endl;
    }

    mynlp = new EndEffectorParametersIdentification();
    mynlp->display_info = display_info;

    app = IpoptApplicationFactory();
    if (mynlp->enable_hessian) {
        app->Options()->SetStringValue("hessian_approximation", "exact");
    }
    else {
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    }
}


void EndEffectorPybindWrapper::set_ipopt_parameters(const double tol,
                                               const double max_wall_time, 
                                               const int print_level,
                                               const int max_iter,
                                               const std::string mu_strategy,
                                               const std::string linear_solver,
                                               const bool gradient_check) {
    app->Options()->SetNumericValue("tol", tol);
    app->Options()->SetNumericValue("max_wall_time", max_wall_time);
    app->Options()->SetIntegerValue("print_level", print_level);
    app->Options()->SetIntegerValue("max_iter", max_iter);
    app->Options()->SetStringValue("mu_strategy", mu_strategy);
    app->Options()->SetStringValue("linear_solver", linear_solver);
    app->Options()->SetStringValue("ma57_automatic_scaling", "yes");

    if (gradient_check) {
        app->Options()->SetStringValue("output_file", "ipopt.out");
        app->Options()->SetStringValue("derivative_test", "first-order");
        app->Options()->SetNumericValue("derivative_test_perturbation", 1e-8);
        app->Options()->SetNumericValue("derivative_test_tol", 1e-5);
        app->Options()->SetNumericValue("point_perturbation_radius", 1);
    }

    set_ipopt_parameters_check = true;
    has_optimized = false;
}

nb::tuple EndEffectorPybindWrapper::optimize() {
    if (!set_ipopt_parameters_check ) {
        throw std::runtime_error("parameters not set properly!");
    }
    // Initialize optimizer
    try {
        mynlp->reset();
        mynlp->set_parameters(model, 
                             trajectory_filename,
                             H,
                             offset);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

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

    has_optimized = mynlp->ifFeasible;

    theta_sol = x_to_theta(mynlp->solution);
    const size_t shape_ptr[] = {theta_sol.size()};
    auto result = nb::ndarray<nb::numpy, const double>(theta_sol.data(),
                                                       1,
                                                       shape_ptr,
                                                       nb::handle());
    return nb::make_tuple(result, mynlp->ifFeasible);
}

double EndEffectorPybindWrapper::analyze_solution() {
    if (!has_optimized) {
        throw std::runtime_error("No optimization has been performed or the optimization is not feasible!");
    }

    Eigen::VectorXd phi_tail = model.inertias[model.nv].toDynamicParameters();
    Eigen::MatrixXd LMI_original = theta_to_LMI(phi_tail);
    Eigen::MatrixXd LMI_solution = theta_to_LMI(theta_sol);
    double distance = BuresWassersteinDistance(LMI_original, LMI_solution);
    std::cout << "Bures-Wasserstein distance: " << distance << std::endl;

    return distance;
}

Eigen::VectorXd EndEffectorPybindWrapper::x_to_theta(Eigen::VectorXd& z) {
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

double EndEffectorPybindWrapper::BuresWassersteinDistance(Eigen::MatrixXd& A, Eigen::MatrixXd& B) {
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

Eigen::MatrixXd EndEffectorPybindWrapper::matrixSquareRoot(Eigen::MatrixXd& matrix) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();

    // Compute the square root of the diagonal matrix of eigenvalues
    Eigen::VectorXd sqrtEigenvalues = eigenvalues.array().sqrt();

    // Reconstruct the square root matrix
    return eigenvectors * sqrtEigenvalues.asDiagonal() * eigenvectors.transpose();
}

Eigen::MatrixXd EndEffectorPybindWrapper::theta_to_LMI(Eigen::VectorXd& theta) {
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
    
}