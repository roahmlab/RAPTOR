#include "KinovaWaitrPybindWrapper.h"

namespace RAPTOR {
namespace Kinova {

KinovaWaitrPybindWrapper::KinovaWaitrPybindWrapper(const std::string urdf_filename) {
    // Define robot model
    pinocchio::urdf::buildModel(urdf_filename, model);

    actual_model_nq = model.nq - 1;

    model.gravity.linear()(2) = GRAVITY;

    qdes = VecX::Zero(actual_model_nq);
    joint_limits_buffer = VecX::Zero(actual_model_nq);
    velocity_limits_buffer = VecX::Zero(actual_model_nq);
    torque_limits_buffer = VecX::Zero(actual_model_nq);

    mynlp = new KinovaWaitrOptimizer();
    app = IpoptApplicationFactory();

    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
}

void KinovaWaitrPybindWrapper::set_obstacles(const nb_2d_double obstacles_inp) {
    if (obstacles_inp.shape(1) != 9) {
        throw std::invalid_argument("Obstacles must have 9 columns, xyz, rpy, size");
    }

    num_obstacles = obstacles_inp.shape(0);

    boxCenters.resize(num_obstacles);
    boxOrientation.resize(num_obstacles);
    boxSize.resize(num_obstacles);

    for (int i = 0; i < num_obstacles; i++) {
        boxCenters[i] << obstacles_inp(i, 0), obstacles_inp(i, 1), obstacles_inp(i, 2);
        boxOrientation[i] << obstacles_inp(i, 3), obstacles_inp(i, 4), obstacles_inp(i, 5);
        boxSize[i] << obstacles_inp(i, 6), obstacles_inp(i, 7), obstacles_inp(i, 8);
    }

    set_obstacles_check = true;
}

void KinovaWaitrPybindWrapper::set_contact_surface_parameters(const double mu_inp,
                                                              const double R_inp,
                                                              const double maxSuctionForce_inp) {
    csp.mu = mu_inp;
    csp.R = R_inp;
    csp.maxSuctionForce = maxSuctionForce_inp;
}

void KinovaWaitrPybindWrapper::set_ipopt_parameters(const double tol,
                                                    const double obj_scaling_factor,
                                                    const double max_wall_time, 
                                                    const int print_level,
                                                    const std::string mu_strategy,
                                                    const std::string linear_solver,
                                                    const bool gradient_check) {
    app->Options()->SetNumericValue("tol", tol);
    app->Options()->SetNumericValue("obj_scaling_factor", obj_scaling_factor);
    app->Options()->SetNumericValue("max_wall_time", max_wall_time);
    app->Options()->SetIntegerValue("print_level", print_level);
    app->Options()->SetStringValue("mu_strategy", mu_strategy);
    app->Options()->SetStringValue("linear_solver", linear_solver);
    app->Options()->SetStringValue("ma57_automatic_scaling", "yes");

    if (gradient_check) {
        app->Options()->SetStringValue("output_file", "ipopt.out");
        app->Options()->SetStringValue("derivative_test", "first-order");
        app->Options()->SetNumericValue("derivative_test_perturbation", 1e-8);
        app->Options()->SetNumericValue("derivative_test_tol", 1e-5);
    }

    set_ipopt_parameters_check = true;
}

void KinovaWaitrPybindWrapper::set_trajectory_parameters(const nb_1d_double q0_inp,
                                                         const nb_1d_double qd0_inp,
                                                         const nb_1d_double qdd0_inp,
                                                         const double duration_inp) {
    if (q0_inp.shape(0) != actual_model_nq || 
        qd0_inp.shape(0) != actual_model_nq || 
        qdd0_inp.shape(0) != actual_model_nq) {
        throw std::invalid_argument("q0, qd0, qdd0 must be of size actual_model_nq");
    }

    T = duration_inp;

    if (T <= 0.0) {
        throw std::invalid_argument("Duration must be positive");
    }

    atp.q0.resize(actual_model_nq);
    atp.q_d0.resize(actual_model_nq);
    atp.q_dd0.resize(actual_model_nq);

    for (int i = 0; i < actual_model_nq; i++) {
        atp.q0(i) = q0_inp(i);
        atp.q_d0(i) = qd0_inp(i);
        atp.q_dd0(i) = qdd0_inp(i);
    }     

    set_trajectory_parameters_check = true;                             
}

void KinovaWaitrPybindWrapper::set_buffer(const nb_1d_double joint_limits_buffer_inp,
                                          const nb_1d_double velocity_limits_buffer_inp,
                                          const nb_1d_double torque_limits_buffer_inp) {
    if (joint_limits_buffer_inp.shape(0) != actual_model_nq || 
        velocity_limits_buffer_inp.shape(0) != actual_model_nq || 
        torque_limits_buffer_inp.shape(0) != actual_model_nq) {
        throw std::invalid_argument("joint_limits_buffer, velocity_limits_buffer, torque_limits_buffer must be of size actual_model_nq");
    }

    for (int i = 0; i < actual_model_nq; i++) {
        joint_limits_buffer(i) = joint_limits_buffer_inp(i);
        velocity_limits_buffer(i) = velocity_limits_buffer_inp(i);
        torque_limits_buffer(i) = torque_limits_buffer_inp(i);
    }                               

    set_buffer_check = true;                                                    
}

void KinovaWaitrPybindWrapper::set_target(const nb_1d_double q_des_inp,
                                          const double tplan_inp) {
    tplan = tplan_inp;

    if (tplan <= 0.0 || tplan > T) {
        throw std::invalid_argument("tplan must be greater than 0.0 or smaller than duration");
    }

    if (q_des_inp.shape(0) != actual_model_nq) {
        throw std::invalid_argument("q_des must be of size actual_model_nq");
    }

    for (int i = 0; i < actual_model_nq; i++) {
        qdes(i) = q_des_inp(i);
    }

    tplan_n = int(tplan / T * N);

    set_target_check = true;
}

nb::tuple KinovaWaitrPybindWrapper::optimize() {
    if (!set_obstacles_check || 
        !set_ipopt_parameters_check || 
        !set_trajectory_parameters_check || 
        !set_buffer_check ||
        !set_target_check) {
        throw std::runtime_error("parameters not set properly!");
    }

    // Define initial guess
    Eigen::VectorXd z(actual_model_nq);
    z.setZero();

    // Initialize Kinova optimizer
    try {
	    mynlp->set_parameters(z,
                              T,
                              N,
                              degree,
                              model,
                              atp,
                              csp,
                              boxCenters,
                              boxOrientation,
                              boxSize,
                              qdes,
                              tplan_n,
                              joint_limits_buffer,
                              velocity_limits_buffer,
                              torque_limits_buffer);
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

    set_ipopt_parameters_check = false;
    set_trajectory_parameters_check = false;
    set_target_check = false;

    const size_t shape_ptr[] = {actual_model_nq};
    auto result = nb::ndarray<nb::numpy, const double>(mynlp->solution.data(),
                                                       1,
                                                       shape_ptr,
                                                       nb::handle());
    return nb::make_tuple(result, mynlp->ifFeasible);
}

}; // namespace Kinova
}; // namespace RAPTOR
