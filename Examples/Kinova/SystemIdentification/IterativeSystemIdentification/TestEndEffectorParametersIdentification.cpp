#include "EndEffectorParametersIdentification.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/rnea.hpp"

using namespace RAPTOR;

int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Error: No arguments provided. please choose the downsampled file" << std::endl;
        return 1; 
    }

    // Initialize model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);

    model.friction.setZero();
    model.damping.setZero();
    model.armature.setZero(); 

    // Initialize parameters
    bool include_offset_input = false;
    Eigen::VectorXd full_parameters = Eigen::VectorXd::Zero((model.nv * (include_offset_input ? 14 : 13)));

    // read inertia parameters
    for (int i = 0; i < model.nv; i++) {
        const int pinocchio_joint_id = i + 1;
        full_parameters.segment<10>(10 * i) =model.inertias[pinocchio_joint_id]
                .toDynamicParameters();
    }

    // read friction parameters
    const std::string solFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/friction_parameters_solution_6.csv";
    Eigen::VectorXd fricition_parameters = Utils::initializeEigenMatrixFromFile(solFile).col(0);

    full_parameters.segment(10 * model.nv, model.nv) = fricition_parameters.head(model.nv);
    full_parameters.segment(11 * model.nv, model.nv)= fricition_parameters.segment(model.nv, model.nv);
    full_parameters.segment(12 * model.nv, model.nv) = fricition_parameters.segment(2 * model.nv, model.nv);
    if (include_offset_input) {
        full_parameters.segment(13 * model.nv, model.nv) = fricition_parameters.tail(model.nv);
    }
 
    // load pos, vel ,acc, tau data from lab
    const std::string posFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/q_downsampled_" 
                                +std::string(argv[1])+ ".csv";
    const std::string velFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/q_d_downsampled_"
                                + std::string(argv[1])+ ".csv";
    const std::string accFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/q_dd_downsampled_"
                                + std::string(argv[1]) +".csv";
    const std::string torqueFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/tau_downsampled_"
                                + std::string(argv[1]) +".csv";


    Eigen::MatrixXd posData = Utils::initializeEigenMatrixFromFile(posFile);
    Eigen::MatrixXd velData = Utils::initializeEigenMatrixFromFile(velFile);
    Eigen::MatrixXd accData = Utils::initializeEigenMatrixFromFile(accFile);
    Eigen::MatrixXd torqueData = Utils::initializeEigenMatrixFromFile(torqueFile);
    const int N = posData.rows();

    // save as sharepoints
    std::shared_ptr<Eigen::VectorXd> full_parametersPtr_ = std::make_shared<Eigen::VectorXd>(full_parameters);
    std::shared_ptr<Eigen::MatrixXd> posDataPtr_ = std::make_shared<Eigen::MatrixXd>(posData.transpose());
    std::shared_ptr<Eigen::MatrixXd> velDataPtr_ = std::make_shared<Eigen::MatrixXd>(velData.transpose());
    std::shared_ptr<Eigen::MatrixXd> accDataPtr_ = std::make_shared<Eigen::MatrixXd>(accData.transpose());
    std::shared_ptr<Eigen::MatrixXd> torqueDataPtr_ = std::make_shared<Eigen::MatrixXd>(torqueData.transpose());

    // Initialize Kinova optimizer
    SmartPtr<EndEffectorParametersIdentification> mynlp = new  EndEffectorParametersIdentification();
    try {
	    mynlp->set_parameters(model,
                              posDataPtr_,
                              velDataPtr_,
                              accDataPtr_,
                              torqueDataPtr_,
                              full_parametersPtr_,
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

        // Write the friction parameters into the file 
        const std::string outputfolder = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/";
        std::ofstream solution(outputfolder + "Endeffector_parameters_solution_" + std::string(argv[1]) + ".csv");

        solution << std::setprecision(16);
        for (int i = 0; i < mynlp->solution.size(); i++) {
            solution << mynlp->solution(i) << std::endl;
        }

        // write the estimate tau into the file  
        const std::string outputfolder1 = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/";
        std::ofstream estimate_tau(outputfolder1 + "Endeffector_estimate_tau_"  + std::string(argv[1]) +".csv");

        // updata the full_parameters from the solution 
        full_parameters.segment(10 * (model.nv -1) , 10) = mynlp ->solution; 

        // inertia tau 
        Eigen::MatrixXd FullObservationMatrix = Eigen::MatrixXd::Zero(N * model.nv, 10 * model.nv);
        for (int i = 0; i < N; i++) {
            const Eigen::VectorXd& q_d = velDataPtr_->col(i);
            const Eigen::VectorXd& q_dd = accDataPtr_->col(i);
            const Eigen::VectorXd& q = posDataPtr_->col(i);

            pinocchio::computeJointTorqueRegressor(
                model, data,
                q, q_d, q_dd);

            FullObservationMatrix.middleRows(i * model.nv, model.nv) = 
                data.jointTorqueRegressor;
        }
        Eigen::VectorXd tau_inertials = FullObservationMatrix * full_parameters.segment(0, 10 * model.nv);

        // friction tau 
        Eigen::VectorXd friction =  full_parameters.segment(10 * model.nv , model.nv);
        Eigen::VectorXd damping = full_parameters.segment(11 * model.nv, model.nv);
        Eigen::VectorXd armature = full_parameters.segment(12 * model.nv, model.nv);
        Eigen::VectorXd offset = Eigen::VectorXd::Zero(model.nv);
        if (include_offset_input) {
            Eigen::VectorXd offset = full_parameters.tail(model.nv);
        }

        for (int i =0 ; i < N; i++) {
            const Eigen::VectorXd& q_d = velDataPtr_->col(i);
            const Eigen::VectorXd& q_dd = accDataPtr_->col(i);
            const Eigen::VectorXd& tau = torqueDataPtr_->col(i);
            const Eigen::VectorXd& q = posDataPtr_->col(i);
            const Eigen::VectorXd& tau_inertial = tau_inertials.segment(i * model.nv, model.nv);
            Eigen::VectorXd tau_estimated = Eigen::VectorXd::Zero(model.nv);

            Eigen::VectorXd total_friction_force = 
                friction.cwiseProduct(q_d.cwiseSign()) +
                damping.cwiseProduct(q_d) +
                armature.cwiseProduct(q_dd) +
                offset;
            tau_estimated = tau_inertial + total_friction_force;
    
            for (int j = 0; j < tau_estimated.size(); ++j) {
                estimate_tau << tau_estimated(j);
                if (j < tau_estimated.size() - 1) {
                    estimate_tau << " ";  
                }
            }
            estimate_tau << std::endl;
        }
    }
    catch (std::exception& e) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    return 0;
}
