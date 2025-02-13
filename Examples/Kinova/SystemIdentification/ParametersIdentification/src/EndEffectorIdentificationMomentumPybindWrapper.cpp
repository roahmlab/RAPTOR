#include "EndEffectorIdentificationMomentumPybindWrapper.h"

namespace RAPTOR {
namespace Kinova {

EndEffectorIdentificationMomentumPybindWrapper::EndEffectorIdentificationMomentumPybindWrapper(const std::string urdf_filename,
                                                                                               const nb_1d_double friction_parameters_input,
                                                                                               const int H_input,
                                                                                               const std::string time_format_string,
                                                                                               const bool display_info):
    H(H_input) {
    // define robot model
    pinocchio::Model model_double;
    pinocchio::urdf::buildModel(urdf_filename, model_double);
    model = model_double.cast<double>();

    // load friction parameters
    if (friction_parameters_input.size() != 3 * model.nv &&
        friction_parameters_input.size() != 4 * model.nv) {
        std::cerr << "Friction_parameters size: " << friction_parameters_input.size() << std::endl;
        throw std::runtime_error("Friction solution input is wrong");
    }

    Eigen::VectorXd friction_parameters = Eigen::Map<Eigen::VectorXd>(friction_parameters_input.data(), friction_parameters_input.size());
    model.friction = friction_parameters.head(model.nv);
    model.damping = friction_parameters.segment(model.nv, model.nv);
    model.armature  = friction_parameters.segment(2 * model.nv, model.nv);
    offset = Eigen::VectorXd::Zero(model.nv);
    if (friction_parameters.size() == 4 * model.nv) {
        offset = friction_parameters.tail(model.nv);
    }

    // time format info
    if (time_format_string == "Second" || 
        time_format_string == "second") {
        time_format = TimeFormat::Second;
    }
    else if (time_format_string == "Millisecond" || 
             time_format_string == "millisecond") {
        time_format = TimeFormat::Millisecond;
    }
    else if (time_format_string == "Microsecond" ||
             time_format_string == "microsecond") {
        time_format = TimeFormat::Microsecond;
    }
    else if (time_format_string == "Nanosecond" ||
             time_format_string == "nanosecond") {
        time_format = TimeFormat::Nanosecond;
    }
    else {
        throw std::runtime_error("Invalid time format string!");
    }

    // sensor noise info
    sensor_noise = SensorNoiseInfo(model.nv);
    sensor_noise.position_error.setZero();
    sensor_noise.velocity_error.setZero();
    sensor_noise.acceleration_error_type = SensorNoiseInfo::SensorNoiseType::Ratio;
    sensor_noise.acceleration_error.setConstant(0.025);

    // set up the Ipopt problem
    mynlp = new EndEffectorParametersIdentificationMomentum();
    mynlp->display_info = display_info;
    mynlp->set_parameters(model,
                          offset);

    app = IpoptApplicationFactory();
}

void EndEffectorIdentificationMomentumPybindWrapper::set_ipopt_parameters(const double tol,
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
}

void EndEffectorIdentificationMomentumPybindWrapper::add_trajectory_file(const std::string filename_input) {
    try {
        auto start = std::chrono::high_resolution_clock::now();

        mynlp->add_trajectory_file(filename_input,
                                   sensor_noise,
                                   H,
                                   time_format,
                                   downsample_rate);

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Regressor computation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds.\n";
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error adding trajectory file! Check previous error message!");
    }
}

nb::tuple EndEffectorIdentificationMomentumPybindWrapper::optimize() {
    if (!set_ipopt_parameters_check) {
        throw std::runtime_error("parameters not set properly!");
    }

    if (mynlp->Aseg.size() == 0) {
        throw std::runtime_error("Trajectory data not added yet!");
    }

    std::cout << "The following trajectory data files are considered:\n";
    for (const auto& filename: mynlp->trajectoryFilenames_) {
        std::cout << "    Trajectory file: " << filename << std::endl;
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
        std::cout << "System identification solve time: " << solve_time << " milliseconds.\n";
    }
    catch (std::exception& e) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    std::cout << "    theta solution: " << mynlp->theta_solution.transpose() << std::endl;
    std::cout << "    theta uncertainties: " << mynlp->theta_uncertainty.transpose() << std::endl;

    const size_t shape_ptr[] = {10};
    auto result1 = nb::ndarray<nb::numpy, const double>(mynlp->theta_solution.data(),
                                                        1,
                                                        shape_ptr,
                                                        nb::handle());

    auto result2 = nb::ndarray<nb::numpy, const double>(mynlp->theta_uncertainty.data(),
                                                        1,
                                                        shape_ptr,
                                                        nb::handle()); 

    return nb::make_tuple(result1, result2);
}

void EndEffectorIdentificationMomentumPybindWrapper::reset() {
    mynlp->reset();
}
    
}; // namespace Kinova
}; // namespace RAPTOR