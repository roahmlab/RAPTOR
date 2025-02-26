#include "KinovaIKSolver.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

using namespace RAPTOR;
using namespace Kinova;
using namespace Ipopt;

int main() {
    // Define robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    
    // Compute forward kinematics at a random configuration first
    ForwardKinematicsSolver fkSolver(&model);

    std::srand(std::time(nullptr));
    Eigen::VectorXd q = Eigen::VectorXd::Random(model.nq);
    fkSolver.compute(0, model.nq, q);
    const Transform desiredTransform = fkSolver.getTransform();

    // Obstacle information (no obstacle in this example)
    int num_obstacles = 0;
    std::vector<Eigen::Vector3d> boxCenters;
    std::vector<Eigen::Vector3d> boxOrientation;
    std::vector<Eigen::Vector3d> boxSize;
    double collision_buffer = 0;

    // Define initial guess
    Eigen::VectorXd z = Eigen::VectorXd::Random(model.nq);

    // Initialize Kinova optimizer
    SmartPtr<KinovaIKSolver> mynlp = new KinovaIKSolver();
    try {
	    mynlp->set_parameters(z,
                              model,
                              desiredTransform,
                              boxCenters,
                              boxOrientation,
                              boxSize);
        mynlp->constr_viol_tol = 1e-8;
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-4);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-6);
	app->Options()->SetNumericValue("max_wall_time", 2.0);
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

    // For gradient checking
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
    catch (std::exception& e) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    std::cout << q.transpose() << std::endl;
    std::cout << mynlp->solution.transpose() << std::endl;

    std::cout << "Desired transform: " << desiredTransform << std::endl;
    fkSolver.compute(0, model.nq, mynlp->solution);
    std::cout << "Solved transform: " << fkSolver.getTransform() << std::endl;

    return 0;
}
