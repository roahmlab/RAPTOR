#include "EndEffectorParametersIdentification.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/rnea.hpp"

using namespace RAPTOR;

int main() {
    // Initialize model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);

    model.friction.setZero();
    model.damping.setZero();
    model.armature.setZero(); 

    // Initialize end-effector parameters
    // Eigen::VectorXd phi = modelPtr_->inertias[modelPtr_->nv] //modelPtr_->nv is pinocchio_joint_id of end_effector 
    //             .toDynamicParameters();
    const std::string solFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/experiment_filter_data/full_parameters_solution.csv";
    Eigen::VectorXd phi = Utils::initializeEigenMatrixFromFile(solFile).col(0);
    std::shared_ptr<Eigen::VectorXd> phiPtr_ = std::make_shared<Eigen::VectorXd>(phi);

    // Initialize data
    const int N = 1322;
    const std::string posFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/experiment_filter_data/q_downsampled.csv";
    const std::string velFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/experiment_filter_data/q_d_downsampled.csv";
    const std::string accFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/experiment_filter_data/q_dd_downsampled.csv";
    const std::string torqueFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/experiment_filter_data/tau_downsampled.csv";


    Eigen::MatrixXd posData = Utils::initializeEigenMatrixFromFile(posFile);
    Eigen::MatrixXd velData = Utils::initializeEigenMatrixFromFile(velFile);
    Eigen::MatrixXd accData = Utils::initializeEigenMatrixFromFile(accFile);
    Eigen::MatrixXd torqueData = Utils::initializeEigenMatrixFromFile(torqueFile);

    std::shared_ptr<Eigen::MatrixXd> posDataPtr_ = std::make_shared<Eigen::MatrixXd>(posData.transpose());
    std::shared_ptr<Eigen::MatrixXd> velDataPtr_ = std::make_shared<Eigen::MatrixXd>(velData.transpose());
    std::shared_ptr<Eigen::MatrixXd> accDataPtr_ = std::make_shared<Eigen::MatrixXd>(accData.transpose());
    std::shared_ptr<Eigen::MatrixXd> torqueDataPtr_ = std::make_shared<Eigen::MatrixXd>(torqueData.transpose());

    // compute nominal torque from data
    for (Index i = 0; i < N; i++) {
        const Eigen::VectorXd& q = posDataPtr_->col(i);
        const Eigen::VectorXd& v = velDataPtr_->col(i);
        const Eigen::VectorXd& a = accDataPtr_->col(i);
        pinocchio::rnea(model, data, q, v, a);
        torqueDataPtr_->col(i) = data.tau + Eigen::VectorXd::Random(model.nv);
    } 

    // Initialize Kinova optimizer
    SmartPtr<EndEffectorParametersIdentification> mynlp = new  EndEffectorParametersIdentification();
    try {
	    mynlp->set_parameters(model,
                              posDataPtr_,
                              velDataPtr_,
                              accDataPtr_,
                              torqueDataPtr_,
                              phiPtr_,
                              true);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-4);
    app->Options()->SetNumericValue("constr_viol_tol", 1e-5);
    mynlp->constr_viol_tol = 1e-5;
    // app->Options()->SetNumericValue("obj_scaling_factor", 1e-3);
	app->Options()->SetNumericValue("max_wall_time", 100.0);
	app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetIntegerValue("max_iter", 3000);
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
    app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "second-order");
    app->Options()->SetStringValue("derivative_test", "first-order");
    app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
    app->Options()->SetNumericValue("derivative_test_tol", 1e-5);

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
        std::cout << "paramer solution: " << mynlp->solution.transpose() << std::endl;
    }
    catch (std::exception& e) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    return 0;
}
