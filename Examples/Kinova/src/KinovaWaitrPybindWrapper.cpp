#include "KinovaWaitrPybindWrapper.h"

namespace IDTO {
namespace Kinova {

KinovaWaitrPybindWrapper::KinovaWaitrPybindWrapper(const std::string urdf_filename) {
    // set openmp number of threads
    int num_threads = 32; // this number is currently hardcoded
    omp_set_num_threads(num_threads);
    
    // Define robot model
    pinocchio::urdf::buildModel(urdf_filename, model);

    actual_model_nq = model.nq - 1;

    model.gravity.linear()(2) = -9.81;

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    jtype.resize(model.nq);
    jtype << 3, 3, 3, 3, 3, 3, 3, 
             3;

    qdes.resize(actual_model_nq);

    mynlp = new KinovaWaitrOptimizer();
    app = IpoptApplicationFactory();

    app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    joint_limits_buffer = VecX::Zero(actual_model_nq);
    velocity_limits_buffer = VecX::Zero(actual_model_nq);
    torque_limits_buffer = VecX::Zero(actual_model_nq);

    // TODO: define contact surface parameters
    // csp.mu = 0.5;
}

void KinovaWaitrPybindWrapper::set_obstacles(const int num_obstacles_inp,
                                             const py::array_t<double> obstacles_inp) {
    num_obstacles = num_obstacles_inp;      

    if (num_obstacles < 0) {
        throw std::invalid_argument("Number of obstacles must be non-negative");
    }                                

    py::buffer_info obstacles_buf = obstacles_inp.request();

    if (obstacles_buf.size != num_obstacles * (3 + 3 + 3)) {
        throw std::invalid_argument("Number of obstacles does not match the input array size");
    }

    const double* obstacles_ptr = (const double*)obstacles_buf.ptr;

    boxCenters.resize(num_obstacles);
    boxOrientation.resize(num_obstacles);
    boxSize.resize(num_obstacles);

    for (int i = 0; i < num_obstacles; i++) {
        boxCenters(i) << obstacles_ptr[i * 9 + 0],
                         obstacles_ptr[i * 9 + 1],
                         obstacles_ptr[i * 9 + 2];
        boxOrientation(i) << obstacles_ptr[i * 9 + 3],
                             obstacles_ptr[i * 9 + 4],
                             obstacles_ptr[i * 9 + 5];
        boxSize(i) << obstacles_ptr[i * 9 + 6],
                      obstacles_ptr[i * 9 + 7],
                      obstacles_ptr[i * 9 + 8];
    }

    set_obstacles_check = true;
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

    if (gradient_check) {
        app->Options()->SetStringValue("output_file", "ipopt.out");
        app->Options()->SetStringValue("derivative_test", "first-order");
        app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
        app->Options()->SetNumericValue("derivative_test_tol", 1e-6);
    }

    set_ipopt_parameters_check = true;
}

void KinovaWaitrPybindWrapper::set_trajectory_parameters(const py::array_t<double> q0_inp,
                                                    const py::array_t<double> qd0_inp,
                                                    const py::array_t<double> qdd0_inp,
                                                    const double duration_inp) {
    py::buffer_info q0_buf = q0_inp.request();
    py::buffer_info qd0_buf = qd0_inp.request();
    py::buffer_info qdd0_buf = qdd0_inp.request();

    T = duration_inp;

    if (q0_buf.size != actual_model_nq) {
        throw std::invalid_argument("q0 must be of size actual_model_nq");
    }

    if (qd0_buf.size != actual_model_nq) {
        throw std::invalid_argument("qd0 must be of size actual_model_nq");
    }
    
    if (qdd0_buf.size != actual_model_nq) {
        throw std::invalid_argument("qdd0 must be of size actual_model_nq");
    }

    if (T <= 0.0) {
        throw std::invalid_argument("Duration must be positive");
    }

    const double *q0_ptr = (const double*)q0_buf.ptr;
    const double *qd0_ptr = (const double*)qd0_buf.ptr;
    const double *qdd0_ptr = (const double*)qdd0_buf.ptr;
    
    atp.q0.resize(actual_model_nq);
    atp.q_d0.resize(actual_model_nq);
    atp.q_dd0.resize(actual_model_nq);

    for (int i = 0; i < actual_model_nq; i++) {
        atp.q0(i) = q0_ptr[i];
        atp.q_d0(i) = qd0_ptr[i];
        atp.q_dd0(i) = qdd0_ptr[i];
    }     

    set_trajectory_parameters_check = true;                             
}

void KinovaWaitrPybindWrapper::set_buffer(const py::array_t<double> joint_limits_buffer_inp,
                                     const py::array_t<double> velocity_limits_buffer_inp,
                                     const py::array_t<double> torque_limits_buffer_inp) {
    py::buffer_info joint_limits_buffer_buf = joint_limits_buffer_inp.request();
    py::buffer_info velocity_limits_buffer_buf = velocity_limits_buffer_inp.request();
    py::buffer_info torque_limits_buffer_buf = torque_limits_buffer_inp.request();

    if (joint_limits_buffer_buf.size != actual_model_nq) {
        throw std::invalid_argument("joint_limits_buffer must be of size actual_model_nq");
    }

    if (velocity_limits_buffer_buf.size != actual_model_nq) {
        throw std::invalid_argument("velocity_limits_buffer must be of size actual_model_nq");
    }

    if (torque_limits_buffer_buf.size != actual_model_nq) {
        throw std::invalid_argument("torque_limits_buffer must be of size actual_model_nq");
    }

    const double* joint_limits_buffer_ptr = (const double*)joint_limits_buffer_buf.ptr;
    const double* velocity_limits_buffer_ptr = (const double*)velocity_limits_buffer_buf.ptr;
    const double* torque_limits_buffer_ptr = (const double*)torque_limits_buffer_buf.ptr;

    for (int i = 0; i < actual_model_nq; i++) {
        joint_limits_buffer(i) = joint_limits_buffer_ptr[i];
        velocity_limits_buffer(i) = velocity_limits_buffer_ptr[i];
        torque_limits_buffer(i) = torque_limits_buffer_ptr[i];
    }                               

    set_buffer_check = true;                                                    
}

void KinovaWaitrPybindWrapper::set_target(const py::array_t<double> q_des_inp,
                                     const double tplan_inp) {
    tplan = tplan_inp;

    if (tplan <= 0.0 || tplan > T) {
        throw std::invalid_argument("tplan must be greater than 0.0 or smaller than duration");
    }

    py::buffer_info q_des_buf = q_des_inp.request();

    if (q_des_buf.size != actual_model_nq) {
        throw std::invalid_argument("q_des must be of size actual_model_nq");
    }

    const double *q_des_ptr = (const double*)q_des_buf.ptr;

    for (int i = 0; i < actual_model_nq; i++) {
        qdes(i) = q_des_ptr[i];
    }
    tplan_n = int(tplan / T * N);

    set_target_check = true;
}

py::tuple KinovaWaitrPybindWrapper::optimize() {
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
                              jtype,
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
    catch (int errorCode) {
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
    catch (int errorCode) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    set_ipopt_parameters_check = false;
    set_trajectory_parameters_check = false;
    set_target_check = false;

    py::array_t<double> result(actual_model_nq, mynlp->solution.data());
    return py::make_tuple(result, mynlp->ifFeasible);
}

}; // namespace Kinova
}; // namespace IDTO
