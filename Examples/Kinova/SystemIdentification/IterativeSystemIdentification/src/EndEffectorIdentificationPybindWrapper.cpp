#include "EndEffectorIdentificationPybindWrapper.h"

namespace RAPTOR {
namespace Kinova {

EndEffectorIdentificationPybindWrapper::EndEffectorIdentificationPybindWrapper(const std::string urdf_filename,
                                                                               const std::string trajectory_filename,
                                                                               const std::string friction_parameters_filename,
                                                                               const int H_input,
                                                                               const bool display_info):
    H(H_input),
    trajectory_filename(trajectory_filename) {
    // define robot model
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
    offset = Eigen::VectorXd::Zero(model.nv);
    if (friction_parameters.size() == 4 * model.nv) {
        offset = friction_parameters.tail(model.nv);
    }

    // sensor noise info
    sensor_noise = SensorNoiseInfo(model.nv);
    sensor_noise.position_error.setZero();
    sensor_noise.velocity_error.setZero();
    sensor_noise.acceleration_error_type = SensorNoiseInfo::SensorNoiseType::Ratio;
    sensor_noise.acceleration_error.setConstant(0.05);

    // set up the Ipopt problem
    mynlp = new EndEffectorParametersIdentification();
    mynlp->display_info = display_info;

    app = IpoptApplicationFactory();
}

void EndEffectorIdentificationPybindWrapper::set_ipopt_parameters(const double tol,
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
        app->Options()->SetStringValue("derivative_test", "second-order");
        app->Options()->SetNumericValue("derivative_test_perturbation", 1e-8);
        app->Options()->SetNumericValue("derivative_test_tol", 1e-5);
        app->Options()->SetNumericValue("point_perturbation_radius", 1);
    }

    set_ipopt_parameters_check = true;
    has_optimized = false;
}

nb::tuple EndEffectorIdentificationPybindWrapper::optimize() {
    if (!set_ipopt_parameters_check ) {
        throw std::runtime_error("parameters not set properly!");
    }
    // Initialize optimizer
    try {
        mynlp->reset();
        mynlp->set_parameters(model, 
                              trajectory_filename,
                              sensor_noise,
                              H,
                              downsample_rate,
                              offset);

        if (mynlp->enable_hessian) {
            app->Options()->SetStringValue("hessian_approximation", "exact");
        }
        else {
            app->Options()->SetStringValue("hessian_approximation", "limited-memory");
        }
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

    VecX theta_sol = mynlp->z_to_theta(mynlp->solution);
    const size_t shape_ptr[] = {theta_sol.size()};
    auto result = nb::ndarray<nb::numpy, const double>(theta_sol.data(),
                                                       1,
                                                       shape_ptr,
                                                       nb::handle());
    return nb::make_tuple(result, mynlp->ifFeasible);
}
    
}; // namespace Kinova
}; // namespace RAPTOR