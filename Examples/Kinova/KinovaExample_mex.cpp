#include "KinovaOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

using namespace IDTO;
using namespace Kinova;
using namespace Ipopt;

using std::cout;
using std::endl;

int main() {
    // Define robot model
    const std::string urdf_filename = "/home/roahmlab/Documents/IDTO/Examples/Kinova/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = -9.81;

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    Eigen::VectorXi jtype(model.nq);
    jtype << 3, 3, 3, 3, 3, 3, 3;

    // Define obstacles
    int num_obstacles = 0;

    // Define trajectories
    ArmourTrajectoryParameters atp;
    atp.q0 = Eigen::VectorXd::Zero(model.nq);
    atp.q_d0 = Eigen::VectorXd::Zero(model.nq);
    atp.q_dd0 = Eigen::VectorXd::Zero(model.nq);

    double T = 1;
    int N = 16;
    const int degree = 5;

    // Define target
    Eigen::VectorXd qdes(model.nq);
    const int tplan_n = N / 2;

    // Define initial guess
    Eigen::VectorXd z(model.nq);
    z.setZero();

    // Define limits buffer
    Eigen::VectorXd joint_limits_buffer(model.nq);
    Eigen::VectorXd torque_limits_buffer(model.nq);

    std::ifstream inputstream("/home/roahmlab/Documents/IDTO/build/oracle_input_buffer.txt");
    if (!inputstream.is_open()) {
        throw std::runtime_error("Error reading input files!");
    }
    for (int i = 0; i < model.nq; i++) {
        inputstream >> atp.q0[i];
    }
    for (int i = 0; i < model.nq; i++) {
        inputstream >> atp.q_d0[i];
    }
    for (int i = 0; i < model.nq; i++) {
        inputstream >> atp.q_dd0[i];
    }
    inputstream >> T;
    for (int i = 0; i < model.nq; i++) {
        inputstream >> z[i];
    }
    for (int i = 0; i < model.nq; i++) {
        inputstream >> qdes[i];
    }
    inputstream >> num_obstacles;

    Eigen::Array<Eigen::Vector3d, 1, Eigen::Dynamic> zonotopeCenters(num_obstacles);
    Eigen::Array<Eigen::MatrixXd, 1, Eigen::Dynamic> zonotopeGenerators(num_obstacles);

    for (int i = 0; i < num_obstacles; i++) {
        for (int j = 0; j < 3; j++) {
            inputstream >> zonotopeCenters(i)(j);
        }
        zonotopeGenerators(i) = Eigen::MatrixXd::Zero(3, MAX_OBSTACLE_GENERATOR_NUM);
        for (int j = 0; j < MAX_OBSTACLE_GENERATOR_NUM; j++) {
            for (int k = 0; k < 3; k++) {
                inputstream >> zonotopeGenerators(i)(k, j);
            }
        }
    }

    for (int i = 0; i < model.nq; i++) {
        inputstream >> joint_limits_buffer[i];
    }
    for (int i = 0; i < model.nq; i++) {
        inputstream >> torque_limits_buffer[i];
    }
    
    inputstream.close();

    // Initialize Kinova optimizer
    SmartPtr<KinovaOptimizer> mynlp = new KinovaOptimizer();
    try {
	    mynlp->set_parameters(z,
                              T,
                              N,
                              degree,
                              model,
                              jtype,
                              atp,
                              zonotopeCenters,
                              zonotopeGenerators,
                              qdes,
                              tplan_n,
                              joint_limits_buffer,
                              torque_limits_buffer);
    }
    catch (int errorCode) {
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-4);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-4);
    // app->Options()->SetNumericValue("obj_scaling_factor", 1e-3);
	app->Options()->SetNumericValue("max_wall_time", 0.2);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetIntegerValue("max_iter", 50);
    app->Options()->SetStringValue("mu_strategy", "monotone");
    app->Options()->SetStringValue("linear_solver", "ma57");
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    // For gradient checking
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-5);

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

    // Print the solution
    std::ofstream solution("solution-kinova.txt");
    if (mynlp->solution.size() == mynlp->numVars && mynlp->ifFeasible) {
        solution << std::setprecision(20);
        for (int i = 0; i < mynlp->numVars; i++) {
            solution << mynlp->solution[i] << std::endl;
        }

        // std::ofstream trajectory("trajectory-kinova.txt");
        // trajectory << std::setprecision(20);
        // for (int i = 0; i < NUM_JOINTS; i++) {
        //     for (int j = 0; j < N; j++) {
        //         trajectory << mynlp->trajPtr_->q(j)(i) << ' ';
        //     }
        //     trajectory << std::endl;
        // }
        // for (int i = 0; i < NUM_JOINTS; i++) {
        //     for (int j = 0; j < N; j++) {
        //         trajectory << mynlp->trajPtr_->q_d(j)(i) << ' ';
        //     }
        //     trajectory << std::endl;
        // }
        // for (int i = 0; i < NUM_JOINTS; i++) {
        //     for (int j = 0; j < N; j++) {
        //         trajectory << mynlp->trajPtr_->q_dd(j)(i) << ' ';
        //     }
        //     trajectory << std::endl;
        // }
        // for (int i = 0; i < NUM_JOINTS; i++) {
        //     for (int j = 0; j < N; j++) {
        //         trajectory << mynlp->idPtr_->tau(j)(i) << ' ';
        //     }
        //     trajectory << std::endl;
        // }
        // trajectory.close();
    }
    else {
        solution << -1 << std::endl;
    }
    solution.close();

    return 0;
}