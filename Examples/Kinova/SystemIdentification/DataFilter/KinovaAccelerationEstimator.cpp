#include "DataFilterOptimizer.h" 

using namespace IDTO;
using namespace Kinova;
using namespace Ipopt;

int main () {
    // set openmp number of threads
    int num_threads = 32; // this number is currently hardcoded
    omp_set_num_threads(num_threads);

    // Define Data
    Eigen::VectorXd tspan = Eigen::VectorXd::LinSpaced(100, 0, 10);
    Eigen::MatrixXd q_data = Eigen::MatrixXd::Zero(7, 100);
    Eigen::MatrixXd q_d_data = Eigen::MatrixXd::Zero(7, 100);

    // Define trajectory parameters
    const int degree = 5;
    const int base_frequency = 10;

    // Define initial guess
    Eigen::VectorXd z = 2 * 1.0 * Eigen::VectorXd::Random((2 * degree + 1) * 7).array() - 1.0;

    // Initialize Kinova optimizer
    SmartPtr<DataFilterOptimizer> mynlp = new DataFilterOptimizer();
    try {
	    mynlp->set_parameters(z,
                              tspan,
                              q_data,
                              q_d_data,
                              degree,
                              base_frequency);
    }
    catch (int errorCode) {
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-6);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-6);
	app->Options()->SetNumericValue("max_wall_time", 10.0);
	app->Options()->SetIntegerValue("print_level", 5);
    // app->Options()->SetIntegerValue("max_iter", 50);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", "ma57");
	app->Options()->SetStringValue("hessian_approximation", "exact");

    // // For gradient checking
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "second-order");
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
    catch (int errorCode) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    return 0;   
}