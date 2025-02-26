#include "ExcitingTrajectoryGenerator.h"

using namespace RAPTOR;
using namespace Kinova;
using namespace Ipopt;

int main(int argc, char* argv[]) {
    // Define robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/gen3_2f85_fixed.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = GRAVITY;
    model.friction.setZero();
    model.damping.setZero();
    model.armature.setZero();

    // Define the indices of identifiable parameters for the Kinova robot, determined solely by its kinematic model
    // Refer to README on how to get this information for different robots
    Eigen::VectorXi independent_param_inds(43);
    independent_param_inds << 9, 11, 12, 14, 19, 18, 17, 15, 21, 22, 24, 29, 28, 27, 25, 31, 32, 34, 39, 38, 37, 35, 41, 42, 44, 49, 48, 47, 45, 51, 52, 54, 59, 58, 57, 55, 61, 62, 64, 69, 68, 67, 65;

    // Define trajectory parameters
    const double T = 10.0;
    const int N = 128;
    const int degree = 3;
    const double base_frequency = 2.0 * M_PI / T;

    // start from a specific static configuration
    Eigen::VectorXd q0(model.nv);
    q0 << -1.16396061,  0.34987045, -3.68094828, -1.75466411,  0.19968747, -1.09177839, -0.19841415;
    Eigen::VectorXd q_d0 = Eigen::VectorXd::Zero(model.nv);

    // Define initial guess
    std::srand(static_cast<unsigned int>(time(0)));
    Eigen::VectorXd z = 0.5 * Eigen::VectorXd::Random((2 * degree + 1) * model.nv);

    // Define obstacles
    std::vector<Eigen::Vector3d> boxCenters = {
        Eigen::Vector3d(0.0, 0.0, 0.18), // floor
        Eigen::Vector3d(0.53, 0.49, 0.56), // back wall
        Eigen::Vector3d(-0.39, -0.84, 0.56), // bar near the control
        Eigen::Vector3d(-0.39, -0.17, 0.56), // bar bewteen 10 and 20 change to wall
        Eigen::Vector3d(0.0, 0.0, 1.12), // ceiling
        Eigen::Vector3d(0.47, -0.09, 1.04) // top camera
    };
    std::vector<Eigen::Vector3d> boxOrientations = {
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0)
    };
    std::vector<Eigen::Vector3d> boxSizes = {
        Eigen::Vector3d(5.0, 5.0, 0.01),
        Eigen::Vector3d(5.0, 0.05, 1.12),
        Eigen::Vector3d(0.05, 0.05, 1.12),
        Eigen::Vector3d(0.05, 1.28, 1.28),
        Eigen::Vector3d(5, 5, 0.05),
        Eigen::Vector3d(0.15, 0.15, 0.15)
    };

    // Define limits buffer
    Eigen::VectorXd joint_limits_buffer(model.nq);
    joint_limits_buffer.setConstant(0.02);
    Eigen::VectorXd velocity_limits_buffer(model.nq);
    velocity_limits_buffer.setConstant(0.05);
    Eigen::VectorXd torque_limits_buffer(model.nq);
    torque_limits_buffer.setConstant(0.5);

    // Initialize Kinova optimizer
    SmartPtr<ExcitingTrajectoryGenerator> mynlp = new ExcitingTrajectoryGenerator();
    try {
	    mynlp->set_parameters(z,
                              T,
                              N,
                              degree,
                              base_frequency,
                              q0,
                              q_d0,
                              model,
                              independent_param_inds,
                              boxCenters,
                              boxOrientations,
                              boxSizes,
                              joint_limits_buffer,
                              velocity_limits_buffer,
                              torque_limits_buffer);
        mynlp->constr_viol_tol = 1e-5;
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-6);
    app->Options()->SetNumericValue("constr_viol_tol", mynlp->constr_viol_tol);
	app->Options()->SetNumericValue("max_wall_time", 60.0);
	app->Options()->SetIntegerValue("print_level", 5);
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
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-3);

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
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    // Re-evaluate the solution at a higher resolution and print the results.
    if (mynlp->ifFeasible) {
        std::shared_ptr<Trajectories> traj = std::make_shared<FixedFrequencyFourierCurves>(T, 
                                                                                           5000, 
                                                                                           model.nv, 
                                                                                           TimeDiscretization::Uniform, 
                                                                                           degree,
                                                                                           base_frequency,
                                                                                           q0,
                                                                                           q_d0);
        std::shared_ptr<RegressorInverseDynamics> rid = std::make_shared<RegressorInverseDynamics>(model, 
                                                                                                   traj,
                                                                                                   true);
        rid->compute(mynlp->solution, false);

        if (argc > 1) {
            const std::string outputfolder = "../Examples/Kinova/SystemIdentification/ExcitingTrajectories/data/T10_d3/";
            std::ofstream solution(outputfolder + "exciting-solution-" + std::string(argv[1]) + ".csv");
            std::ofstream trajectory(outputfolder + "exciting-trajectory-" + std::string(argv[1]) + ".csv");

            solution << std::setprecision(16);
            for (int i = 0; i < mynlp->solution.size(); i++) {
                solution << mynlp->solution(i) << std::endl;
            }
            for (int i = 0; i < 7; ++i){
                solution << q0(i) << std::endl;
            }
            for (int i = 0; i < 7; ++i){
                solution << q_d0(i) << std::endl;
            }
            solution << base_frequency << std::endl;

            trajectory << std::setprecision(16);
            for (int i = 0; i < traj->N; i++) {
                trajectory << traj->tspan(i) << ' ';
                trajectory << traj->q(i).transpose() << ' ';
                trajectory << traj->q_d(i).transpose() << ' ';
                trajectory << rid->tau(i).transpose() << std::endl;
            }
        }
        else {
            std::ofstream solution("exciting-solution.csv");
            std::ofstream trajectory("exciting-trajectory.csv");

            solution << std::setprecision(16);
            for (int i = 0; i < mynlp->solution.size(); i++) {
                solution << mynlp->solution(i) << std::endl;
            }
            for (int i = 0; i < 7; ++i){
                solution << q0(i) << std::endl;
            }
            for (int i = 0; i < 7; ++i){
                solution << q_d0(i) << std::endl;
            }
            solution << base_frequency << std::endl;

            for (int i = 0; i < traj->N; i++) {
                trajectory << traj->tspan(i) << ' ';
                trajectory << traj->q(i).transpose() << ' ';
                trajectory << traj->q_d(i).transpose() << ' ';
                trajectory << rid->tau(i).transpose() << std::endl;
            }
        }
    }

    return 0;
}