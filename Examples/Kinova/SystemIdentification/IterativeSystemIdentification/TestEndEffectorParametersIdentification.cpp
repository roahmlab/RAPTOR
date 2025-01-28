#include "EndEffectorParametersIdentification.h"

#include "BezierCurves.h"

using namespace RAPTOR;

int main(int argc, char* argv[]) {
    // Load the robot model
    const std::string urdf_filename  = "../Robots/kinova-gen3/gen3_2f85_fixed.urdf";

    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);

    model.friction.setConstant(1.0);
    model.damping.setConstant(5.0);
    model.armature.setConstant(0.0);
    Eigen::VectorXd offset = Eigen::VectorXd::Zero(model.nv);

    // Generate a random trajectory data
    const double T = 10.0; // duration of the trajectory
    const int N = 1000; // number of samples in the data

        // to compute the acceleration and then the torque
    std::shared_ptr<BezierCurves> test_trajectory = std::make_shared<BezierCurves>(
        T, N, model.nv, TimeDiscretization::Uniform, 3);
    test_trajectory->compute(Eigen::VectorXd::Random(test_trajectory->varLength), false);

        // to store the position, velocity and torque data for system identification
    std::shared_ptr<TrajectoryData> trajectory_data = std::make_shared<TrajectoryData>(
        T, N, model.nv);

        // the matrix that stores the acceleration data
    Eigen::MatrixXd acceleration(N, model.nv);

    for (int i = 0; i < N; i++) {
        acceleration.row(i) = test_trajectory->q_dd(i);

        pinocchio::rnea(model, data, test_trajectory->q(i), test_trajectory->q_d(i), test_trajectory->q_dd(i));

        trajectory_data->q(i) = test_trajectory->q(i);
        trajectory_data->q_d(i) = test_trajectory->q_d(i);
        trajectory_data->q_dd(i) = data.tau + 
                                   model.friction.cwiseProduct(test_trajectory->q_d(i).cwiseSign()) +
                                   model.damping.cwiseProduct(test_trajectory->q_d(i)) +
                                   offset;
    }
  
    // Initialize the Ipopt problem
    SmartPtr<EndEffectorParametersIdentification> mynlp = new EndEffectorParametersIdentification();

    double setup_time = 0;
    try {
        auto start = std::chrono::high_resolution_clock::now();
	    mynlp->set_parameters(model,
                              offset);
        mynlp->add_trajectory_file(trajectory_data,
                                  acceleration);
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
	app->Options()->SetNumericValue("max_wall_time", 1.0);
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
