#include "KinovaLongerHorizonPybindWrapper.h"

namespace RAPTOR {
namespace Kinova {

KinovaLongerHorizonPybindWrapper::KinovaLongerHorizonPybindWrapper(const std::string urdf_filename,
                                                                   const bool display_info) {
    // Define robot model
    pinocchio::Model model_double;
    pinocchio::urdf::buildModel(urdf_filename, model_double);
    model = model_double.cast<double>();
    
    model.gravity.linear()(2) = GRAVITY;

    mynlp = new KinovaLongerHorizonOptimizer();
    mynlp->display_info = display_info;

    app = IpoptApplicationFactory();

    joint_limits_buffer = VecX::Zero(model.nv);
    velocity_limits_buffer = VecX::Zero(model.nv);
    torque_limits_buffer = VecX::Zero(model.nv);
}

void KinovaLongerHorizonPybindWrapper::set_obstacles(const nb_2d_double obstacles_inp,
                                                     const double collision_buffer_inp) {
    if (obstacles_inp.shape(1) != 9) {
        throw std::invalid_argument("Obstacles must have 9 columns, xyz, rpy, size");
    }

    num_obstacles = obstacles_inp.shape(0);
    collision_buffer = collision_buffer_inp;

    boxCenters.resize(num_obstacles);
    boxOrientation.resize(num_obstacles);
    boxSize.resize(num_obstacles);

    for (int i = 0; i < num_obstacles; i++) {
        boxCenters[i] << obstacles_inp(i, 0), obstacles_inp(i, 1), obstacles_inp(i, 2);
        boxOrientation[i] << obstacles_inp(i, 3), obstacles_inp(i, 4), obstacles_inp(i, 5);
        boxSize[i] << obstacles_inp(i, 6), obstacles_inp(i, 7), obstacles_inp(i, 8);
    }

    set_obstacles_check = true;
    has_optimized = false;
}

void KinovaLongerHorizonPybindWrapper::set_ipopt_parameters(const double tol,
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
        app->Options()->SetStringValue("derivative_test", "second-order");
        app->Options()->SetNumericValue("derivative_test_perturbation", 1e-8);
        app->Options()->SetNumericValue("derivative_test_tol", 1e-4);
    }

    set_ipopt_parameters_check = true;
    has_optimized = false;
}

void KinovaLongerHorizonPybindWrapper::set_trajectory_parameters(const nb_1d_double q0_inp,
                                                                 const nb_1d_double qT_inp,
                                                                 const double duration_inp,
                                                                 const int degree_inp) {
    if (q0_inp.shape(0) != model.nv || 
        qT_inp.shape(0) != model.nv) {
        throw std::invalid_argument("q0, qT must be of size model.nv");
    }

    T = duration_inp;
    degree = degree_inp;

    if (T <= 0.0) {
        throw std::invalid_argument("Duration must be positive");
    }

    if (degree < 1) {
        throw std::invalid_argument("Degree must be at least 1");
    }

    q0.resize(model.nv);
    qT.resize(model.nv);

    for (int i = 0; i < model.nv; i++) {
        q0(i) = q0_inp(i);
        qT(i) = qT_inp(i);
    }     

    set_trajectory_parameters_check = true;        
    has_optimized = false;                     
}

void KinovaLongerHorizonPybindWrapper::set_buffer(const nb_1d_double joint_limits_buffer_inp,
                                                  const nb_1d_double velocity_limits_buffer_inp,
                                                  const nb_1d_double torque_limits_buffer_inp) {
    if (joint_limits_buffer_inp.shape(0) != model.nv || 
        velocity_limits_buffer_inp.shape(0) != model.nv || 
        torque_limits_buffer_inp.shape(0) != model.nv) {
        throw std::invalid_argument("joint_limits_buffer, velocity_limits_buffer, torque_limits_buffer must be of size model.nv");
    }

    for (int i = 0; i < model.nv; i++) {
        joint_limits_buffer(i) = joint_limits_buffer_inp(i);
        velocity_limits_buffer(i) = velocity_limits_buffer_inp(i);
        torque_limits_buffer(i) = torque_limits_buffer_inp(i);
    }                               

    set_buffer_check = true;    
    has_optimized = false;                                                
}

nb::tuple KinovaLongerHorizonPybindWrapper::optimize() {
    if (!set_obstacles_check || 
        !set_ipopt_parameters_check || 
        !set_trajectory_parameters_check || 
        !set_buffer_check) {
        throw std::runtime_error("parameters not set properly!");
    }

    // Define initial guess:
    // a series of Bezier curves with degree + 1 control points, straight from q0 to qT
    VecX z0 = VecX::Zero(model.nv * degree * 3);
    VecX qdiff = (qT - q0) / (degree + 1);
    for (int i = 0; i < degree; i++) {
        z0.segment(i * model.nv * 3, model.nv) = q0 + qdiff * (i + 1);
    }

    N = 20 * (degree + 1); // number of time steps to check constraints

    // Initialize Kinova optimizer
    try {
        mynlp->reset();
	    mynlp->set_parameters(z0,
                              T,
                              N,
                              degree,
                              model,
                              q0,
                              qT,
                              boxCenters,
                              boxOrientation,
                              boxSize,
                              joint_limits_buffer,
                              velocity_limits_buffer,
                              torque_limits_buffer,
                              true,
                              collision_buffer);
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

    set_trajectory_parameters_check = false;
    has_optimized = mynlp->ifFeasible;

    const size_t shape_ptr[] = {mynlp->solution.size()};
    auto result = nb::ndarray<nb::numpy, const double>(mynlp->solution.data(),
                                                       1,
                                                       shape_ptr,
                                                       nb::handle());
    return nb::make_tuple(result, mynlp->ifFeasible);
}

nb::tuple KinovaLongerHorizonPybindWrapper::analyze_solution() {
    if (!has_optimized) {
        std::cerr << "Warning: No optimization has been performed or the optimization is not feasible!" << std::endl;
    }

    // re-evaluate the solution on a finer time discretization
    const int N_simulate = T * 24; // replay for 24 Hz
    
    SmartPtr<KinovaLongerHorizonOptimizer> testnlp = new KinovaLongerHorizonOptimizer();
    try {
        testnlp->display_info = false;
        testnlp->set_parameters(mynlp->solution,
                                T,
                                N_simulate,
                                degree,
                                model,
                                q0,
                                qT,
                                boxCenters,
                                boxOrientation,
                                boxSize,
                                joint_limits_buffer,
                                velocity_limits_buffer,
                                torque_limits_buffer,
                                true,
                                collision_buffer);
        Index n, m, nnz_jac_g, nnz_h_lag;
        TNLP::IndexStyleEnum index_style;
        testnlp->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
        Number ztry[testnlp->numVars], x_l[testnlp->numVars], x_u[testnlp->numVars];
        Number g[testnlp->numCons], g_lb[testnlp->numCons], g_ub[testnlp->numCons];
        for (int i = 0; i < testnlp->numVars; i++) {
            ztry[i] = mynlp->solution(i);
        }
        testnlp->get_bounds_info(testnlp->numVars, x_l, x_u, testnlp->numCons, g_lb, g_ub);
        testnlp->eval_g(testnlp->numVars, ztry, false, testnlp->numCons, g);
        testnlp->summarize_constraints(testnlp->numCons, g, false);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error evaluating the solution on a finer time discretization! Check previous error message!");
    }

    trajInfo.resize(N_simulate, 4 * NUM_JOINTS + 1);
    for (int i = 0; i < N_simulate; i++) {
        for (int j = 0; j < NUM_JOINTS; j++) {
            trajInfo(i, j) = testnlp->trajPtr_->q(i)(j);
            trajInfo(i, j + NUM_JOINTS) = testnlp->trajPtr_->q_d(i)(j);
            trajInfo(i, j + NUM_JOINTS*2) = testnlp->trajPtr_->q_dd(i)(j);
            trajInfo(i, j + NUM_JOINTS*3) = testnlp->idPtr_->tau(i)(j);
        }

        trajInfo(i, 4 * NUM_JOINTS) = testnlp->trajPtr_->tspan(i);
    }

    const KinovaCustomizedConstraints* customizedConstraintsPtr = nullptr;
    for (int i = 0; i < testnlp->constraintsNameVec_.size(); i++) {
        if (testnlp->constraintsNameVec_[i] == "obstacle avoidance constraints") {
            customizedConstraintsPtr = dynamic_cast<KinovaCustomizedConstraints*>(testnlp->constraintsPtrVec_[i].get());
        }
    }
    if (customizedConstraintsPtr == nullptr) {
        throw std::runtime_error("Error retrieving customized constraints!");
    }

    spheres_x.resize(customizedConstraintsPtr->num_spheres, N_simulate);
    spheres_y.resize(customizedConstraintsPtr->num_spheres, N_simulate);
    spheres_z.resize(customizedConstraintsPtr->num_spheres, N_simulate);
    for (int i = 0; i < N_simulate; i++) {
        for (int j = 0; j < customizedConstraintsPtr->num_spheres; j++) {
            const Vec3& sphere_center = customizedConstraintsPtr->sphere_centers_copy(j, i);
            spheres_x(j, i) = sphere_center(0);
            spheres_y(j, i) = sphere_center(1);
            spheres_z(j, i) = sphere_center(2);
        }
    }

    const size_t shape_ptr1[] = {N_simulate, 4 * NUM_JOINTS + 1};
    auto traj = nb::ndarray<nb::numpy, double, nb::shape<2, -1>>(
        trajInfo.data(),
        2,
        shape_ptr1,
        nb::handle()
    );

    const size_t shape_ptr2[] = {customizedConstraintsPtr->num_spheres, N_simulate};
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

    return nb::make_tuple(traj, sphere_xs, sphere_ys, sphere_zs);
}

}; // namespace Kinova
}; // namespace RAPTOR
