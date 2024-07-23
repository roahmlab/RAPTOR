#include "ConditionNumberOptimizer.h"

using namespace RAPTOR;
using namespace Kinova;
using namespace Ipopt;

int main(int argc, char* argv[]) {
    // set openmp number of threads
    int num_threads = 32; // this number is currently hardcoded
    omp_set_num_threads(num_threads);

    const std::string regroupMatrixFileName = "../Examples/Kinova/SystemIdentification/RegroupMatrix.csv";
    // const std::string regroupMatrixFileName = "../Examples/Kinova/SystemIdentification/RegroupMatrix_withoutmotordynamics.csv";
    
    // Define robot model
    const std::string urdf_filename = "../Examples/Kinova/ArmourUnsafe/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    model.gravity.linear()(2) = -9.81;
    model.friction.setZero();
    // model.damping.setZero();
    // model.rotorInertia.setZero();

    // manually define the joint axis of rotation
    // 1 for Rx, 2 for Ry, 3 for Rz
    // 4 for Px, 5 for Py, 6 for Pz
    // not sure how to extract this from a pinocchio model so define outside here.
    Eigen::VectorXi jtype(model.nq);
    jtype << 3, 3, 3, 3, 3, 3, 3;

    // Define trajectory parameters
    const double T = 10.0;
    const int N = 50;
    const int degree = 5;
    const double base_frequency = 2.0 * M_PI / T;

    // Define initial guess
    std::srand(static_cast<unsigned int>(time(0)));
    Eigen::VectorXd z = 2 * 0.2 * Eigen::VectorXd::Random((2 * degree + 3) * model.nv).array() - 0.1;
    z.block((2 * degree + 1) * model.nv, 0, model.nv, 1) = 
        2 * 1.0 * Eigen::VectorXd::Random(model.nv).array() - 1.0;
    z.block((2 * degree + 1) * model.nv + model.nv, 0, model.nv, 1) = 
        2 * 0.5 * Eigen::VectorXd::Random(model.nv).array() - 0.5;

    // Define limits buffer
    Eigen::VectorXd joint_limits_buffer(model.nq);
    joint_limits_buffer.setConstant(0.02);
    Eigen::VectorXd velocity_limits_buffer(model.nq);
    velocity_limits_buffer.setConstant(0.05);
    Eigen::VectorXd torque_limits_buffer(model.nq);
    torque_limits_buffer.setConstant(0.5);

    // Initialize Kinova optimizer
    SmartPtr<ConditionNumberOptimizer> mynlp = new ConditionNumberOptimizer();
    try {
	    mynlp->set_parameters(z,
                              T,
                              N,
                              degree,
                              base_frequency,
                              model,
                              jtype,
                              regroupMatrixFileName,
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

    // For gradient checking
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

        const Eigen::VectorXd initial_position_part = mynlp->solution.segment((2 * degree + 1) * model.nv, model.nv);
        mynlp->solution.segment((2 * degree + 1) * model.nv, model.nv) =
            Utils::wrapToPi(initial_position_part);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    // Re-evaluate the solution at a higher resolution and print the results.
    if (mynlp->ifFeasible) {
        std::shared_ptr<Trajectories> traj = std::make_shared<FixedFrequencyFourierCurves>(T, 
                                                                                           2000, 
                                                                                           model.nv, 
                                                                                           TimeDiscretization::Uniform, 
                                                                                           degree,
                                                                                           base_frequency);
        traj->compute(mynlp->solution, false);

        if (argc > 1) {
            const std::string outputfolder = "../Examples/Kinova/SystemIdentification/ExcitingTrajectories/data/T10_d5_slower/";
            std::ofstream solution(outputfolder + "exciting-solution-" + std::string(argv[1]) + ".csv");
            std::ofstream position(outputfolder + "exciting-position-" + std::string(argv[1]) + ".csv");
            std::ofstream velocity(outputfolder + "exciting-velocity-" + std::string(argv[1]) + ".csv");
            std::ofstream acceleration(outputfolder + "exciting-acceleration-" + std::string(argv[1]) + ".csv");

            solution << std::setprecision(16);
            position << std::setprecision(16);
            velocity << std::setprecision(16);
            acceleration << std::setprecision(16);

            for (int i = 0; i < traj->N; i++) {
                position << traj->q(i).transpose() << std::endl;
                velocity << traj->q_d(i).transpose() << std::endl;
                acceleration << traj->q_dd(i).transpose() << std::endl;
            }

            for (int i = 0; i < mynlp->solution.size(); i++) {
                solution << mynlp->solution(i) << std::endl;
            }
            solution << base_frequency << std::endl;
        }
        else {
            std::ofstream solution("exciting-solution.csv");
            std::ofstream position("exciting-position.csv");
            std::ofstream velocity("exciting-velocity.csv");
            std::ofstream acceleration("exciting-acceleration.csv");

            for (int i = 0; i < traj->N; i++) {
                position << traj->q(i).transpose() << std::endl;
                velocity << traj->q_d(i).transpose() << std::endl;
                acceleration << traj->q_dd(i).transpose() << std::endl;
            }

            for (int i = 0; i < mynlp->solution.size(); i++) {
                solution << mynlp->solution(i) << std::endl;
            }
        }
    }

    return 0;
}