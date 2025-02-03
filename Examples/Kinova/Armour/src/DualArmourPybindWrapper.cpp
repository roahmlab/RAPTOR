#include "DualArmourPybindWrapper.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

DualArmourPybindWrapper::DualArmourPybindWrapper(const std::string urdf_filename1,
                                                 const std::string config_filename1,
                                                 const std::string urdf_filename2,
                                                 const std::string config_filename2,
                                                 const bool display_info) {
    robotInfoPtr1_ = std::make_shared<RobotInfo>(urdf_filename1, config_filename1);
    robotInfoPtr2_ = std::make_shared<RobotInfo>(urdf_filename2, config_filename2);

    q_des = VecX::Zero(robotInfoPtr1_->num_motors + robotInfoPtr2_->num_motors);

    mynlp = new DualArmourOptimizer();
    mynlp->display_info = display_info;

    app = IpoptApplicationFactory();
    if (mynlp->enable_hessian) {
        app->Options()->SetStringValue("hessian_approximation", "exact");
    }
    else {
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    }
}

void DualArmourPybindWrapper::set_endeffector_inertial_parameters_robot1(const double object_mass,
                                                                         const nb_1d_double object_com,
                                                                         const nb_1d_double object_inertia) {
    if (object_inertia.shape(0) != 9 || 
        object_com.shape(0) != 3) {
        throw std::invalid_argument("Object inertia must have 9 elements and object com must have 3 elements");
    }
    
    Mat3 object_inertia_mat;
    object_inertia_mat << object_inertia(0), object_inertia(1), object_inertia(2),
                          object_inertia(3), object_inertia(4), object_inertia(5),
                          object_inertia(6), object_inertia(7), object_inertia(8);

    robotInfoPtr1_->model.inertias[robotInfoPtr1_->model.nbodies - 1] = pinocchio::Inertia(
        object_mass, 
        Vec3(object_com(0), 
             object_com(1), 
             object_com(2)),
        object_inertia_mat);
    
    if (trajPtr1_ != nullptr) {
        if (dynPtr1_ == nullptr) {
            dynPtr1_ = std::make_shared<PZDynamics>(robotInfoPtr1_, trajPtr1_);
        }
        else {
            dynPtr1_->reset_robot_info(robotInfoPtr1_);
        }
    }
    else {
        throw std::runtime_error("Trajectory parameters for robot 1 not set yet!");
    }
}

void DualArmourPybindWrapper::set_endeffector_inertial_parameters_robot2(const double object_mass,
                                                                         const nb_1d_double object_com,
                                                                         const nb_1d_double object_inertia) {
    if (object_inertia.shape(0) != 9 || 
        object_com.shape(0) != 3) {
        throw std::invalid_argument("Object inertia must have 9 elements and object com must have 3 elements");
    }
    
    Mat3 object_inertia_mat;
    object_inertia_mat << object_inertia(0), object_inertia(1), object_inertia(2),
                          object_inertia(3), object_inertia(4), object_inertia(5),
                          object_inertia(6), object_inertia(7), object_inertia(8);

    robotInfoPtr2_->model.inertias[robotInfoPtr2_->model.nbodies - 1] = pinocchio::Inertia(
        object_mass, 
        Vec3(object_com(0), 
             object_com(1), 
             object_com(2)),
        object_inertia_mat);
    
    if (trajPtr2_ != nullptr) {
        if (dynPtr2_ == nullptr) {
            dynPtr2_ = std::make_shared<PZDynamics>(robotInfoPtr2_, trajPtr2_);
        }
        else {
            dynPtr2_->reset_robot_info(robotInfoPtr2_);
        }
    }
    else {
        throw std::runtime_error("Trajectory parameters for robot 2 not set yet!");
    }
}

void DualArmourPybindWrapper::set_obstacles(const nb_2d_double obstacles_inp) {
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

void DualArmourPybindWrapper::set_ipopt_parameters(const double tol,
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

void DualArmourPybindWrapper::set_trajectory_parameters(const nb_1d_double q0_inp,
                                                        const nb_1d_double q_d0_inp,
                                                        const nb_1d_double q_dd0_inp,
                                                        const nb_1d_double k_center_inp,
                                                        const nb_1d_double k_range_inp,
                                                        const double duration_inp,
                                                        const nb_1d_double q_des_inp,
                                                        const double t_plan_inp) {
    if (q0_inp.shape(0) != robotInfoPtr1_->num_motors + robotInfoPtr2_->num_motors || 
        q_d0_inp.shape(0) != robotInfoPtr1_->num_motors + robotInfoPtr2_->num_motors || 
        q_dd0_inp.shape(0) != robotInfoPtr1_->num_motors + robotInfoPtr2_->num_motors) {
        std::cerr << q0_inp.shape(0) << " " << q_d0_inp.shape(0) << " " << q_dd0_inp.shape(0) << std::endl;
        std::cerr << robotInfoPtr1_->num_motors << " " << robotInfoPtr2_->num_motors << std::endl;
        throw std::invalid_argument("q0, q_d0, q_dd0 must be of size total motors");
    }

    if (k_center_inp.shape(0) != robotInfoPtr1_->num_motors + robotInfoPtr2_->num_motors || 
        k_range_inp.shape(0) != robotInfoPtr1_->num_motors + robotInfoPtr2_->num_motors) {
        std::cerr << k_center_inp.shape(0) << " " << k_range_inp.shape(0) << std::endl;
        std::cerr << robotInfoPtr1_->num_motors + robotInfoPtr2_->num_motors << std::endl;
        throw std::invalid_argument("k_center, k_range must be of size total motors");
    }

    if (duration_inp <= 0.0) {
        throw std::invalid_argument("Duration must be positive");
    }

    if (q_des_inp.shape(0) != robotInfoPtr1_->num_motors + robotInfoPtr2_->num_motors) {
        std::cerr << q_des_inp.shape(0) << std::endl;
        std::cerr << robotInfoPtr1_->num_motors + robotInfoPtr2_->num_motors << std::endl;
        throw std::invalid_argument("q_des must be of size total motors");
    }

    if (t_plan_inp <= 0.0 || t_plan_inp > 1.0) {
        throw std::invalid_argument("t_plan must be positive and less than 1.0");
    }

    q0_robot1.resize(robotInfoPtr1_->num_motors);
    q_d0_robot1.resize(robotInfoPtr1_->num_motors);
    q_dd0_robot1.resize(robotInfoPtr1_->num_motors);
    for (size_t i = 0; i < robotInfoPtr1_->num_motors; i++) {
        q0_robot1(i) = q0_inp(i);
        q_d0_robot1(i) = q_d0_inp(i);
        q_dd0_robot1(i) = q_dd0_inp(i);
    }
    q0_robot2.resize(robotInfoPtr2_->num_motors);
    q_d0_robot2.resize(robotInfoPtr2_->num_motors);
    q_dd0_robot2.resize(robotInfoPtr2_->num_motors);
    for (size_t i = 0; i < robotInfoPtr2_->num_motors; i++) {
        q0_robot2(i) = q0_inp(robotInfoPtr1_->num_motors + i);
        q_d0_robot2(i) = q_d0_inp(robotInfoPtr1_->num_motors + i);
        q_dd0_robot2(i) = q_dd0_inp(robotInfoPtr1_->num_motors + i);
    }

    k_center_robot1.resize(robotInfoPtr1_->num_motors);
    k_range_robot1.resize(robotInfoPtr1_->num_motors);
    for (size_t i = 0; i < robotInfoPtr1_->num_motors; i++) {
        k_center_robot1(i) = k_center_inp(i);
        k_range_robot1(i) = k_range_inp(i);
    }
    k_center_robot2.resize(robotInfoPtr2_->num_motors);
    k_range_robot2.resize(robotInfoPtr2_->num_motors);
    for (size_t i = 0; i < robotInfoPtr2_->num_motors; i++) {
        k_center_robot2(i) = k_center_inp(robotInfoPtr1_->num_motors + i);
        k_range_robot2(i) = k_range_inp(robotInfoPtr1_->num_motors + i);
    }

    duration = duration_inp;  

    q_des.resize(robotInfoPtr1_->num_motors + robotInfoPtr2_->num_motors);
    for (size_t i = 0; i < robotInfoPtr1_->num_motors; i++) {
        q_des(i) = q_des_inp(i);
    }
    for (size_t i = 0; i < robotInfoPtr2_->num_motors; i++) {
        q_des(robotInfoPtr1_->num_motors + i) = q_des_inp(i);
    }

    t_plan = t_plan_inp;

    trajPtr1_ = std::make_shared<BezierCurveInterval>(
        q0_robot1, q_d0_robot1, q_dd0_robot1, 
        k_center_robot1, k_range_robot1, 
        duration, 
        robotInfoPtr1_);

    if (dynPtr1_ == nullptr) {
        dynPtr1_ = std::make_shared<PZDynamics>(robotInfoPtr1_, trajPtr1_);
    }
    else {
        dynPtr1_->reset_trajectory(trajPtr1_);
    }

    trajPtr2_ = std::make_shared<BezierCurveInterval>(
        q0_robot2, q_d0_robot2, q_dd0_robot2, 
        k_center_robot2, k_range_robot2, 
        duration, 
        robotInfoPtr2_);

    if (dynPtr2_ == nullptr) {
        dynPtr2_ = std::make_shared<PZDynamics>(robotInfoPtr2_, trajPtr2_);
    }
    else {
        dynPtr2_->reset_trajectory(trajPtr2_);
    }

    set_trajectory_parameters_check = true;        
    has_optimized = false;                     
}

nb::tuple DualArmourPybindWrapper::optimize() {
    if (!set_obstacles_check || 
        !set_ipopt_parameters_check || 
        !set_trajectory_parameters_check) {
        throw std::runtime_error("parameters not set properly!");
    }

    // generate Joint Trajectory Reachable Sets
    auto start1 = std::chrono::high_resolution_clock::now();
    dynPtr1_->compute();
    dynPtr2_->compute();
    auto end1 = std::chrono::high_resolution_clock::now();
    std::cout << "    DualArmourPybindWrapper: Time taken to generate reachable sets: " 
              << std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count() / 1000.0
              << " ms" << std::endl;

    // Initialize Kinova optimizer
    try {
        mynlp->reset();
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
        solve_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
        std::cout << "Total solve time: " << solve_time << " milliseconds.\n";
    }
    catch (std::exception& e) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    set_trajectory_parameters_check = false;
    has_optimized = mynlp->ifFeasible;

    const size_t shape_ptr1[] = {robotInfoPtr1_->num_motors};
    auto result1 = nb::ndarray<nb::numpy, const double>(mynlp->solution.data(),
                                                        1,
                                                        shape_ptr1,
                                                        nb::handle());

    const size_t shape_ptr2[] = {robotInfoPtr2_->num_motors};
    auto result2 = nb::ndarray<nb::numpy, const double>(mynlp->solution.data() + robotInfoPtr1_->num_motors,
                                                        1,
                                                        shape_ptr2,
                                                        nb::handle());

    return nb::make_tuple(result1, result2, mynlp->ifFeasible);
}

nb::tuple DualArmourPybindWrapper::analyze_solution() {
    if (!has_optimized) {
        // throw std::runtime_error("No optimization has been performed or the optimization is not feasible!");
        std::cerr << "DualArmourPybindWrapper::Warning: No optimization has been performed or the optimization is not feasible!" << std::endl;
    }
    
    // recover trajectory information
    VecX q(robotInfoPtr1_->num_motors);
    VecX qd(robotInfoPtr1_->num_motors);
    VecX qdd(robotInfoPtr1_->num_motors);

    trajInfo_robot1.resize(3 * robotInfoPtr1_->num_motors + 1, trajPtr1_->num_time_steps);
    for (size_t i = 0; i < trajPtr1_->num_time_steps; i++) {
        const double t = duration * (i + 0.5) / trajPtr1_->num_time_steps;
        trajPtr1_->computeTrajectories(mynlp->solution.head(robotInfoPtr1_->num_motors), t, q, qd, qdd);
        trajInfo_robot1.col(i).head(robotInfoPtr1_->num_motors) = q;
        trajInfo_robot1.col(i).segment(robotInfoPtr1_->num_motors, robotInfoPtr1_->num_motors) = qd;
        trajInfo_robot1.col(i).segment(robotInfoPtr1_->num_motors * 2, robotInfoPtr1_->num_motors) = qdd;
        trajInfo_robot1(robotInfoPtr1_->num_motors * 3, i) = t;
    }

    trajInfo_robot2.resize(3 * robotInfoPtr2_->num_motors + 1, trajPtr2_->num_time_steps);
    for (size_t i = 0; i < trajPtr2_->num_time_steps; i++) {
        const double t = duration * (i + 0.5) / trajPtr2_->num_time_steps;
        trajPtr2_->computeTrajectories(mynlp->solution.tail(robotInfoPtr2_->num_motors), t, q, qd, qdd);
        trajInfo_robot2.col(i).head(robotInfoPtr2_->num_motors) = q;
        trajInfo_robot2.col(i).segment(robotInfoPtr2_->num_motors, robotInfoPtr2_->num_motors) = qd;
        trajInfo_robot2.col(i).segment(robotInfoPtr2_->num_motors * 2, robotInfoPtr2_->num_motors) = qdd;
        trajInfo_robot2(robotInfoPtr2_->num_motors * 3, i) = t;
    }

    // // recover sphere occupancy information for collision checking
    // spheres_x.resize(robotInfoPtr_->num_spheres, trajPtr_->num_time_steps);
    // spheres_y.resize(robotInfoPtr_->num_spheres, trajPtr_->num_time_steps);
    // spheres_z.resize(robotInfoPtr_->num_spheres, trajPtr_->num_time_steps);
    // spheres_radius.resize(robotInfoPtr_->num_spheres, trajPtr_->num_time_steps);
    // for (size_t i = 0; i < trajPtr_->num_time_steps; i++) {
    //     for (size_t j = 0; j < robotInfoPtr_->num_spheres; j++) {
    //         const std::string sphere_name = "collision-" + std::to_string(j);
    //         const pinocchio::FrameIndex frame_id = 
    //             robotInfoPtr_->model.getFrameId(sphere_name);

    //         const auto& PZsphere = dynPtr_->data_sparses[i].oMf[frame_id].translation();

    //         const Interval x_res = PZsphere(0).slice(mynlp->solution);
    //         const Interval y_res = PZsphere(1).slice(mynlp->solution);
    //         const Interval z_res = PZsphere(2).slice(mynlp->solution);
            
    //         spheres_x(j, i) = getCenter(x_res);
    //         spheres_y(j, i) = getCenter(y_res);
    //         spheres_z(j, i) = getCenter(z_res);
    //         spheres_radius(j, i) = dynPtr_->sphere_radii(j, i);
    //     }
    // }

    // // recover torque reachable set information
    // torque_center.resize(robotInfoPtr_->num_motors, trajPtr_->num_time_steps);
    // for (size_t i = 0; i < trajPtr_->num_time_steps; i++) {
    //     for (size_t j = 0; j < robotInfoPtr_->num_motors; j++) {
    //         const PZSparse& PZtorque = dynPtr_->data_sparses[i].tau(j);
    //         const Interval res = PZtorque.slice(mynlp->solution);
    //         torque_center(j, i) = getCenter(res);
    //     }
    // }

    // torque_radius = dynPtr_->torque_radii;

    // // recover contact constraints (reachable sets) here
    // if (robotInfoPtr_->num_joints - robotInfoPtr_->num_motors == 1) {
    //     separation_force_center.resize(trajPtr_->num_time_steps);
    //     separation_force_radius.resize(trajPtr_->num_time_steps);
    //     friction_cone_center.resize(FRICTION_CONE_LINEARIZED_SIZE, trajPtr_->num_time_steps);
    //     friction_cone_radius.resize(FRICTION_CONE_LINEARIZED_SIZE, trajPtr_->num_time_steps);
    //     zmp_center.resize(ZMP_LINEARIZED_SIZE, trajPtr_->num_time_steps);
    //     zmp_radius.resize(ZMP_LINEARIZED_SIZE, trajPtr_->num_time_steps);

    //     for (size_t i = 0; i < trajPtr_->num_time_steps; i++) {
    //         // (1) support force larger than 0
    //         const Interval supportForceRange = 
    //             dynPtr_->data_sparses[i].f[dynPtr_->model_sparses[i].nv].linear()(2)
    //                 .slice(mynlp->solution);
    //         separation_force_center(i) = getCenter(supportForceRange);
    //         separation_force_radius(i) = getRadius(supportForceRange);

    //         // (2) friction cone constraints
    //         for (size_t j = 0; j < FRICTION_CONE_LINEARIZED_SIZE; j++) {
    //             const Interval frictionConstraintRange = 
    //                 dynPtr_->friction_PZs(j, i).slice(mynlp->solution);
    //             friction_cone_center(j, i) = getCenter(frictionConstraintRange);
    //             friction_cone_radius(j, i) = getRadius(frictionConstraintRange);
    //         }      

    //         // (3) ZMP constraints
    //         for (size_t j = 0; j < ZMP_LINEARIZED_SIZE; j++) {
    //             const Interval zmpConstraintRange = 
    //                 dynPtr_->zmp_PZs(j, i).slice(mynlp->solution);
    //             zmp_center(j, i) = getCenter(zmpConstraintRange);
    //             zmp_radius(j, i) = getRadius(zmpConstraintRange);
    //         }
    //     }
    // }

    // export to outputs
    const size_t shape_ptr1[] = {3 * robotInfoPtr1_->num_motors + 1, trajPtr1_->num_time_steps};
    auto traj_robot1 = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
        trajInfo_robot1.data(),
        2,
        shape_ptr1,
        nb::handle()
    );

    const size_t shape_ptr2[] = {3 * robotInfoPtr2_->num_motors + 1, trajPtr2_->num_time_steps};
    auto traj_robot2 = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
        trajInfo_robot2.data(),
        2,
        shape_ptr2,
        nb::handle()
    );

    return nb::make_tuple(traj_robot1, traj_robot2);

    // const size_t shape_ptr1[] = {3 * robotInfoPtr_->num_motors + 1, trajPtr_->num_time_steps};
    // auto traj = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
    //     trajInfo_robot1.data(),
    //     2,
    //     shape_ptr1,
    //     nb::handle()
    // );

    // const size_t shape_ptr2[] = {robotInfoPtr_->num_spheres, trajPtr_->num_time_steps};
    // auto sphere_xs = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
    //     spheres_x.data(),
    //     2,
    //     shape_ptr2,
    //     nb::handle()
    // );
    // auto sphere_ys = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
    //     spheres_y.data(),
    //     2,
    //     shape_ptr2,
    //     nb::handle()
    // );
    // auto sphere_zs = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
    //     spheres_z.data(),
    //     2,
    //     shape_ptr2,
    //     nb::handle()
    // );
    // auto sphere_radii = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
    //     spheres_radius.data(),
    //     2,
    //     shape_ptr2,
    //     nb::handle()
    // );

    // const size_t shape_ptr3[] = {robotInfoPtr_->num_motors, trajPtr_->num_time_steps};
    // auto torque_centers = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
    //     torque_center.data(),
    //     2,
    //     shape_ptr3,
    //     nb::handle()
    // );
    // auto torque_radii = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
    //     torque_radius.data(),
    //     2,
    //     shape_ptr3,
    //     nb::handle()
    // );

    // const size_t shape_ptr4[] = {trajPtr_->num_time_steps};
    // auto separation_force_centers = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
    //     separation_force_center.data(),
    //     1,
    //     shape_ptr4,
    //     nb::handle()
    // );
    // auto separation_force_radii = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
    //     separation_force_radius.data(),
    //     1,
    //     shape_ptr4,
    //     nb::handle()
    // );

    // const size_t shape_ptr5[] = {FRICTION_CONE_LINEARIZED_SIZE, trajPtr_->num_time_steps};
    // auto friction_cone_centers = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
    //     friction_cone_center.data(),
    //     2,
    //     shape_ptr5,
    //     nb::handle()
    // );
    // auto friction_cone_radii = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
    //     friction_cone_radius.data(),
    //     2,
    //     shape_ptr5,
    //     nb::handle()
    // );

    // const size_t shape_ptr6[] = {ZMP_LINEARIZED_SIZE, trajPtr_->num_time_steps};
    // auto zmp_centers = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
    //     zmp_center.data(),
    //     2,
    //     shape_ptr6,
    //     nb::handle()
    // );
    // auto zmp_radii = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
    //     zmp_radius.data(),
    //     2,
    //     shape_ptr6,
    //     nb::handle()
    // );

    // if (robotInfoPtr_->num_joints - robotInfoPtr_->num_motors == 1) {
    //     return nb::make_tuple(traj, sphere_xs, sphere_ys, sphere_zs, sphere_radii, torque_centers, torque_radii, 
    //                           separation_force_centers, separation_force_radii, 
    //                           friction_cone_centers, friction_cone_radii, 
    //                           zmp_centers, zmp_radii);
    // }
    // return nb::make_tuple(traj, sphere_xs, sphere_ys, sphere_zs, sphere_radii, torque_centers, torque_radii);
}

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR
