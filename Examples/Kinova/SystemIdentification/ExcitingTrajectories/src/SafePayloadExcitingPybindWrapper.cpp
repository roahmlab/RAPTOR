#include "SafePayloadExcitingPybindWrapper.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

SafePayloadExcitingPybindWrapper::SafePayloadExcitingPybindWrapper(const std::string urdf_filename,
                                                                   const std::string config_filename,
                                                                   const bool display_info) {
    robotInfoPtr_ = std::make_shared<RobotInfo>(urdf_filename, config_filename);

    mynlp = new SafePayloadExcitingTrajectoryGenerator();
    mynlp->display_info = display_info;

    app = IpoptApplicationFactory();
    if (mynlp->enable_hessian) {
        app->Options()->SetStringValue("hessian_approximation", "exact");
    }
    else {
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    }
}

void SafePayloadExcitingPybindWrapper::set_obstacles(const nb_2d_double obstacles_inp) {
    if (obstacles_inp.shape(1) != 9) {
        throw std::invalid_argument("Obstacles must have 9 columns, xyz, rpy, size");
    }

    num_obstacles = obstacles_inp.shape(0);

    boxCenters.resize(num_obstacles);
    boxOrientation.resize(num_obstacles);
    boxSize.resize(num_obstacles);

    for (size_t i = 0; i < num_obstacles; i++) {
        boxCenters[i] << obstacles_inp(i, 0), obstacles_inp(i, 1), obstacles_inp(i, 2);
        boxOrientation[i] << obstacles_inp(i, 3), obstacles_inp(i, 4), obstacles_inp(i, 5);
        boxSize[i] << obstacles_inp(i, 6), obstacles_inp(i, 7), obstacles_inp(i, 8);
    }

    set_obstacles_check = true;
    has_optimized = false;
}

void SafePayloadExcitingPybindWrapper::set_ipopt_parameters(const double tol,
                                                            const double constr_viol_tol,
                                                            const double obj_scaling_factor,
                                                            const double max_wall_time, 
                                                            const int print_level,
                                                            const std::string mu_strategy,
                                                            const std::string linear_solver,
                                                            const bool gradient_check) {
    app->Options()->SetNumericValue("tol", tol);
    app->Options()->SetNumericValue("constr_viol_tol", constr_viol_tol);
    mynlp->constr_viol_tol = constr_viol_tol;
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
        app->Options()->SetNumericValue("derivative_test_tol", 1e-4);
    }

    set_ipopt_parameters_check = true;
    has_optimized = false;
}

void SafePayloadExcitingPybindWrapper::set_trajectory_parameters(const nb_1d_double q0_inp,
                                                                 const nb_1d_double q_d0_inp,
                                                                 const nb_1d_double q_dd0_inp,
                                                                 const nb_1d_double k_center_inp,
                                                                 const nb_1d_double k_range_inp,
                                                                 const double duration_inp) {
    if (q0_inp.shape(0) != robotInfoPtr_->num_motors || 
        q_d0_inp.shape(0) != robotInfoPtr_->num_motors || 
        q_dd0_inp.shape(0) != robotInfoPtr_->num_motors) {
        std::cerr << q0_inp.shape(0) << " " << q_d0_inp.shape(0) << " " << q_dd0_inp.shape(0) << std::endl;
        std::cerr << robotInfoPtr_->num_motors << std::endl;
        throw std::invalid_argument("q0, q_d0, q_dd0 must be of size robotInfoPtr_->num_motors");
    }

    if (k_center_inp.shape(0) != robotInfoPtr_->num_motors || 
        k_range_inp.shape(0) != robotInfoPtr_->num_motors) {
        std::cerr << k_center_inp.shape(0) << " " << k_range_inp.shape(0) << std::endl;
        std::cerr << robotInfoPtr_->num_motors << std::endl;
        throw std::invalid_argument("k_center, k_range must be of size robotInfoPtr_->num_motors");
    }

    if (duration_inp <= 0.0) {
        throw std::invalid_argument("Duration must be positive");
    }

    q0.resize(robotInfoPtr_->num_motors);
    q_d0.resize(robotInfoPtr_->num_motors);
    q_dd0.resize(robotInfoPtr_->num_motors);
    for (size_t i = 0; i < robotInfoPtr_->num_motors; i++) {
        q0(i) = q0_inp(i);
        q_d0(i) = q_d0_inp(i);
        q_dd0(i) = q_dd0_inp(i);
    }   

    k_center.resize(robotInfoPtr_->num_motors);
    k_range.resize(robotInfoPtr_->num_motors);
    for (size_t i = 0; i < robotInfoPtr_->num_motors; i++) {
        k_center(i) = k_center_inp(i);
        k_range(i) = k_range_inp(i);
    }

    duration = duration_inp;

    trajPtr_ = std::make_shared<BezierCurveInterval>(
        q0, q_d0, q_dd0, 
        k_center, k_range, 
        duration, 
        robotInfoPtr_);

    if (dynPtr_ == nullptr) {
        dynPtr_ = std::make_shared<PZDynamics>(robotInfoPtr_, trajPtr_);
    }
    else {
        dynPtr_->reset_trajectory(trajPtr_);
    }

    set_trajectory_parameters_check = true;        
    has_optimized = false;                     
}

void SafePayloadExcitingPybindWrapper::set_endeffector_inertial_parameters(const nb_1d_double inertial_parameters,
                                                                           const nb_1d_double inertial_parameters_lb,
                                                                           const nb_1d_double inertial_parameters_ub) {
    if (inertial_parameters.shape(0) != 10 || 
        inertial_parameters_lb.shape(0) != 10 ||
        inertial_parameters_ub.shape(0) != 10) {
        throw std::invalid_argument("Inertial parameters must be of size 10");
    }

    Vec10 inertial_parameters_vec;
    Vec10 inertial_parameters_lb_vec;
    Vec10 inertial_parameters_ub_vec;
    for (size_t i = 0; i < 10; i++) {
        inertial_parameters_vec(i) = inertial_parameters(i);
        inertial_parameters_lb_vec(i) = inertial_parameters_lb(i);
        inertial_parameters_ub_vec(i) = inertial_parameters_ub(i);
    }

    robotInfoPtr_->change_endeffector_inertial_parameters(inertial_parameters_vec, 
                                                          inertial_parameters_lb_vec,
                                                          inertial_parameters_ub_vec);

    if (trajPtr_ != nullptr) {
        if (dynPtr_ == nullptr) {
            dynPtr_ = std::make_shared<PZDynamics>(robotInfoPtr_, trajPtr_);
        }
        else {
            dynPtr_->reset_robot_info(robotInfoPtr_);
        }
    }
}

nb::tuple SafePayloadExcitingPybindWrapper::optimize() {
    if (!set_obstacles_check || 
        !set_ipopt_parameters_check || 
        !set_trajectory_parameters_check) {
        throw std::runtime_error("parameters not set properly!");
    }

    // generate Joint Trajectory Reachable Sets
    auto start1 = std::chrono::high_resolution_clock::now();
    dynPtr_->compute();
    auto end1 = std::chrono::high_resolution_clock::now();
    std::cout << "    SafePayloadExcitingPybindWrapper: Time taken to generate reachable sets: " 
              << std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count() / 1000.0
              << " ms" << std::endl;

    double t_plan = 1.0 * duration;

    // Initialize Kinova optimizer
    try {
        mynlp->reset();
	    mynlp->set_parameters(robotInfoPtr_, 
                              trajPtr_, 
                              dynPtr_, 
                              boxCenters,
                              boxOrientation,
                              boxSize);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if ( status != Solve_Succeeded ) {
		throw std::runtime_error("Error during initialization of optimization!");
    }

    // Run ipopt to solve the optimization problem
    double solve_time = 0;
    try {
        auto start = std::chrono::high_resolution_clock::now();

        // Ask Ipopt to solve the problem
        status = app->OptimizeTNLP(mynlp);

        auto end = std::chrono::high_resolution_clock::now();
        solve_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
        std::cout << "Total solve time: " << solve_time << " milliseconds.\n";
    }
    catch (std::exception& e) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    if (status == Invalid_Problem_Definition) {
        // reachable set size too large
        // usually is torque reachable sets, print the size out for debugging
        Utils::writeEigenMatrixToFile(dynPtr_->torque_radii, "torque_radii_exciting.txt");
        throw std::runtime_error("Invalid problem definition!");
    }

    set_trajectory_parameters_check = false;
    has_optimized = mynlp->ifFeasible;

    const size_t shape_ptr[] = {mynlp->solution.size()};
    auto result = nb::ndarray<nb::numpy, const double>(mynlp->solution.data(),
                                                       1,
                                                       shape_ptr,
                                                       nb::handle());
    return nb::make_tuple(result, mynlp->ifFeasible);
}

nb::tuple SafePayloadExcitingPybindWrapper::analyze_solution() {
    if (!has_optimized) {
        throw std::runtime_error("No optimization has been performed or the optimization is not feasible!");
    }
    
    // recover trajectory information
    trajInfo.resize(3 * robotInfoPtr_->num_motors + 1, trajPtr_->num_time_steps);
    for (size_t i = 0; i < trajPtr_->num_time_steps; i++) {
        VecX q(robotInfoPtr_->num_motors);
        VecX qd(robotInfoPtr_->num_motors);
        VecX qdd(robotInfoPtr_->num_motors);
        const double t = duration * (i + 0.5) / trajPtr_->num_time_steps;
        trajPtr_->computeTrajectories(mynlp->solution, t, q, qd, qdd);
        trajInfo.col(i).head(robotInfoPtr_->num_motors) = q;
        trajInfo.col(i).segment(robotInfoPtr_->num_motors, robotInfoPtr_->num_motors) = qd;
        trajInfo.col(i).segment(robotInfoPtr_->num_motors * 2, robotInfoPtr_->num_motors) = qdd;
        trajInfo(robotInfoPtr_->num_motors * 3, i) = t;
    }

    // recover sphere occupancy information for collision checking
    spheres_x.resize(robotInfoPtr_->num_spheres, trajPtr_->num_time_steps);
    spheres_y.resize(robotInfoPtr_->num_spheres, trajPtr_->num_time_steps);
    spheres_z.resize(robotInfoPtr_->num_spheres, trajPtr_->num_time_steps);
    spheres_radius.resize(robotInfoPtr_->num_spheres, trajPtr_->num_time_steps);
    for (size_t i = 0; i < trajPtr_->num_time_steps; i++) {
        for (size_t j = 0; j < robotInfoPtr_->num_spheres; j++) {
            const std::string sphere_name = "collision-" + std::to_string(j);
            const pinocchio::FrameIndex frame_id = 
                robotInfoPtr_->model.getFrameId(sphere_name);

            const auto& PZsphere = dynPtr_->data_sparses[i].oMf[frame_id].translation();

            const Interval x_res = PZsphere(0).slice(mynlp->solution);
            const Interval y_res = PZsphere(1).slice(mynlp->solution);
            const Interval z_res = PZsphere(2).slice(mynlp->solution);
            
            spheres_x(j, i) = getCenter(x_res);
            spheres_y(j, i) = getCenter(y_res);
            spheres_z(j, i) = getCenter(z_res);
            spheres_radius(j, i) = dynPtr_->sphere_radii(j, i);
        }
    }

    // recover torque reachable set information
    torque_center.resize(robotInfoPtr_->num_motors, trajPtr_->num_time_steps);
    for (size_t i = 0; i < trajPtr_->num_time_steps; i++) {
        for (size_t j = 0; j < robotInfoPtr_->num_motors; j++) {
            const PZSparse& PZtorque = dynPtr_->data_sparses[i].tau(j);
            const Interval res = PZtorque.slice(mynlp->solution);
            torque_center(j, i) = getCenter(res);
        }
    }

    torque_radius = dynPtr_->torque_radii;

    // export to outputs
    const size_t shape_ptr1[] = {3 * robotInfoPtr_->num_motors + 1, trajPtr_->num_time_steps};
    auto traj = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
        trajInfo.data(),
        2,
        shape_ptr1,
        nb::handle()
    );

    const size_t shape_ptr2[] = {robotInfoPtr_->num_spheres, trajPtr_->num_time_steps};
    auto sphere_xs = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
        spheres_x.data(),
        2,
        shape_ptr2,
        nb::handle()
    );
    auto sphere_ys = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
        spheres_y.data(),
        2,
        shape_ptr2,
        nb::handle()
    );
    auto sphere_zs = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
        spheres_z.data(),
        2,
        shape_ptr2,
        nb::handle()
    );
    auto sphere_radii = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
        spheres_radius.data(),
        2,
        shape_ptr2,
        nb::handle()
    );

    const size_t shape_ptr3[] = {robotInfoPtr_->num_motors, trajPtr_->num_time_steps};
    auto torque_centers = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
        torque_center.data(),
        2,
        shape_ptr3,
        nb::handle()
    );
    auto torque_radii = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
        torque_radius.data(),
        2,
        shape_ptr3,
        nb::handle()
    );

    return nb::make_tuple(traj, sphere_xs, sphere_ys, sphere_zs, sphere_radii, torque_centers, torque_radii);
}

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR
