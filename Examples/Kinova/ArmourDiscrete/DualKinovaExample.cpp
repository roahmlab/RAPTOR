#include "DualKinovaOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

using namespace RAPTOR;
using namespace Kinova;
using namespace Ipopt;

int main() {
    // Define robot model
    const std::string urdf_filename1 = "../Robots/kinova-gen3/gen3_2f85_fixed.urdf";
    pinocchio::Model model1;
    pinocchio::urdf::buildModel(urdf_filename1, model1);
    model1.gravity.linear()(2) = GRAVITY;

    const std::string urdf_filename2 = "../Robots/kinova-gen3/gen3_2f85_fixed_the_other_side.urdf";
    pinocchio::Model model2;
    pinocchio::urdf::buildModel(urdf_filename2, model2);
    model2.gravity.linear()(2) = GRAVITY;

    // Define obstacles
    const int num_obstacles = 5;
    std::vector<Eigen::Vector3d> boxCenters;
    std::vector<Eigen::Vector3d> boxOrientation;
    std::vector<Eigen::Vector3d> boxSize;

    boxCenters.resize(num_obstacles);
    boxOrientation.resize(num_obstacles);
    boxSize.resize(num_obstacles);

    for (int i = 0; i < num_obstacles; i++) {
        boxCenters[i] << 10, 0, 0.3 * i;
        boxOrientation[i] << 0, 0, 0;
        boxSize[i] << 0.1, 0.1, 0.1;
    }

    // Define trajectories
    Eigen::VectorXd q0 = Eigen::VectorXd::Zero(model1.nv + model2.nv);
    Eigen::VectorXd qT = Eigen::VectorXd::Zero(model1.nv + model2.nv);

    const double T = 10.0;
    const int N = 20;
    // const int degree = 4;

    // Define initial guess
    Eigen::VectorXd z = Eigen::VectorXd::Zero(7 * 4 * 3 * 2);

    // Define limits buffer
    Eigen::VectorXd joint_limits_buffer(model1.nq);
    joint_limits_buffer.setConstant(0.0);
    Eigen::VectorXd velocity_limits_buffer(model1.nq);
    velocity_limits_buffer.setConstant(0.0);
    Eigen::VectorXd torque_limits_buffer(model1.nq);
    torque_limits_buffer.setConstant(0.0);

    // Initialize Kinova optimizer
    SmartPtr<DualKinovaOptimizer> mynlp = new DualKinovaOptimizer();
    try {
	    mynlp->set_parameters(z,
                              T,
                              N,
                              model1,
                              model2,
                              q0,
                              qT,
                              boxCenters,
                              boxOrientation,
                              boxSize,
                              joint_limits_buffer,
                              velocity_limits_buffer,
                              torque_limits_buffer);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-4);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-4);
    // app->Options()->SetNumericValue("obj_scaling_factor", 1e-3);
	app->Options()->SetNumericValue("max_wall_time", 20.0);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetIntegerValue("max_iter", 500);
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
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-4);

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

    // // Print the solution
    // if (mynlp->solution.size() == mynlp->numVars) {
    //     std::ofstream solution("solution-kinova.txt");
    //     solution << std::setprecision(20);
    //     for (int i = 0; i < mynlp->numVars; i++) {
    //         solution << mynlp->solution[i] << std::endl;
    //     }
    //     solution.close();

    //     std::ofstream trajectory("trajectory-kinova.txt");
    //     trajectory << std::setprecision(20);
    //     for (int i = 0; i < NUM_JOINTS; i++) {
    //         for (int j = 0; j < N; j++) {
    //             trajectory << mynlp->trajPtr_->q(j)(i) << ' ';
    //         }
    //         trajectory << std::endl;
    //     }
    //     for (int i = 0; i < NUM_JOINTS; i++) {
    //         for (int j = 0; j < N; j++) {
    //             trajectory << mynlp->trajPtr_->q_d(j)(i) << ' ';
    //         }
    //         trajectory << std::endl;
    //     }
    //     for (int i = 0; i < NUM_JOINTS; i++) {
    //         for (int j = 0; j < N; j++) {
    //             trajectory << mynlp->trajPtr_->q_dd(j)(i) << ' ';
    //         }
    //         trajectory << std::endl;
    //     }
    //     for (int i = 0; i < NUM_JOINTS; i++) {
    //         for (int j = 0; j < N; j++) {
    //             trajectory << mynlp->idPtr_->tau(j)(i) << ' ';
    //         }
    //         trajectory << std::endl;
    //     }
    //     trajectory.close();
    // }

    return 0;
}
