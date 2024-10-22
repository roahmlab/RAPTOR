#include "BaseParametersIdentification.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/rnea.hpp"

using namespace RAPTOR;

int main(int argc, char* argv[]) {
    // Initialize model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);
    pinocchio::Data data1(model);

    model.friction.setZero();
    model.damping.setZero();
    model.armature.setZero(); 

    // Initialize Regroup matrix
    auto qrSolverPtr_ = std::make_shared<QRDecompositionSolver>(model); 
    qrSolverPtr_->generateRandomObservation(4000);
    qrSolverPtr_->computeRegroupMatrix();



    // Initialize data
    bool include_offset_input = true ;
    const std::string posFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/q_downsampled_5.csv";
    const std::string velFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/q_d_downsampled_5.csv";
    const std::string accFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/q_dd_downsampled_5.csv";
    const std::string torqueFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/tau_downsampled_5.csv";


    Eigen::MatrixXd posData = Utils::initializeEigenMatrixFromFile(posFile);
    Eigen::MatrixXd velData = Utils::initializeEigenMatrixFromFile(velFile);
    Eigen::MatrixXd accData = Utils::initializeEigenMatrixFromFile(accFile);
    Eigen::MatrixXd torqueData = Utils::initializeEigenMatrixFromFile(torqueFile);
    torqueData = -torqueData;    

    int N = posData.rows();

    std::shared_ptr<Eigen::MatrixXd> posDataPtr_ = std::make_shared<Eigen::MatrixXd>(posData.transpose());
    std::shared_ptr<Eigen::MatrixXd> velDataPtr_ = std::make_shared<Eigen::MatrixXd>(velData.transpose());
    std::shared_ptr<Eigen::MatrixXd> accDataPtr_ = std::make_shared<Eigen::MatrixXd>(accData.transpose());
    std::shared_ptr<Eigen::MatrixXd> torqueDataPtr_ = std::make_shared<Eigen::MatrixXd>(torqueData.transpose());

    // std::shared_ptr<Eigen::MatrixXd> posDataPtr_ = std::make_shared<Eigen::MatrixXd>();
    // std::shared_ptr<Eigen::MatrixXd> velDataPtr_ = std::make_shared<Eigen::MatrixXd>();
    // std::shared_ptr<Eigen::MatrixXd> accDataPtr_ = std::make_shared<Eigen::MatrixXd>();
    // std::shared_ptr<Eigen::MatrixXd> torqueDataPtr_ = std::make_shared<Eigen::MatrixXd>();

    // posDataPtr_->resize(model.nv, N);
    // velDataPtr_->resize(model.nv, N);
    // accDataPtr_->resize(model.nv, N);
    // torqueDataPtr_->resize(model.nv, N);

    // posDataPtr_->setRandom();
    // velDataPtr_->setRandom();
    // accDataPtr_->setRandom();
    
    // compute nominal torque from data
    // for (Index i = 0; i < N; i++) {
    //     const Eigen::VectorXd& q = posDataPtr_->col(i);
    //     const Eigen::VectorXd& v = velDataPtr_->col(i);
    //     const Eigen::VectorXd& a = accDataPtr_->col(i);
    //     pinocchio::rnea(model, data, q, v, a);
    //     torqueDataPtr_->col(i) = data.tau + Eigen::VectorXd::Random(model.nv);
    // } 

    // Initialize Kinova optimizer
    SmartPtr<BaseParametersIdentification> mynlp = new BaseParametersIdentification();
    try {
	    mynlp->set_parameters(model,
                              posDataPtr_,
                              velDataPtr_,
                              accDataPtr_,
                              torqueDataPtr_,
                              qrSolverPtr_,
                              include_offset_input);
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
	app->Options()->SetNumericValue("max_wall_time", 120.0);
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
    // app->Options()->SetStringValue("output_file", "ipopt.out");
    // app->Options()->SetStringValue("derivative_test", "second-order");
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
       
        Eigen::VectorXd beta = mynlp->solution.head(qrSolverPtr_->dim_id);
        mynlp->solution.segment(0, 10 * model.nv) = qrSolverPtr_->Ginv * mynlp->solution.segment(0, 10 * model.nv);
        std::cout << "paramer solution: " << (mynlp->solution.head(70)).transpose() << std::endl;
        std::cout <<"from urdf paramer" << qrSolverPtr_->phi.transpose() << std::endl;
        // std::cout <<"different" << qrSolverPtr_->phi.transpose() -mynlp->solution.transpose().segment(0,70) << std::endl;

        const auto& phi_diff = qrSolverPtr_->phi.transpose() - mynlp->solution.transpose().segment(0, 70);
        std::cout <<"from urdf paramer - solution"   << std::endl;
        for (int i = 0; i < phi_diff.size(); ++i) {
            std::cout << phi_diff(i) << " ";

            if ((i + 1) % 10 == 0) {
                std::cout << std::endl;
            }
        }

        const std::string outputfolder = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/";
        std::ofstream solution(outputfolder + "full_parameters_solution_" + std::string(argv[1]) + ".csv");

        solution << std::setprecision(16);
        for (int i = 0; i < mynlp->solution.size(); i++) {
            solution << mynlp->solution(i) << std::endl;
        }


        // estimate torque
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

        Eigen::MatrixXd RegroupedObservationMatrix = 
            FullObservationMatrix * qrSolverPtr_->RegroupMatrix;

        Eigen::VectorXd tau_inertials = RegroupedObservationMatrix * beta;

        Eigen::VectorXd friction =  mynlp->solution.segment(10 * model.nv , model.nv);
        Eigen::VectorXd damping = mynlp->solution.segment(11 * model.nv, model.nv);
        Eigen::VectorXd armature = mynlp->solution.segment(12 * model.nv, model.nv);
        Eigen::VectorXd offset = Eigen::VectorXd::Zero(model.nv);
        if (include_offset_input) {
            Eigen::VectorXd offset = mynlp->solution.tail(model.nv);
        }


        const std::string outputfolder1 = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/";
        std::ofstream estimate_tau(outputfolder1 + "full_parameters_estimate_tau_" + std::string(argv[1]) + ".csv");


        for (Index i =0 ; i < N; i++) {
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

        // // Extract the segment from the solution vector (no transpose involved)
        // Eigen::VectorXd phi_solution_segment = mynlp->solution.segment(0, 70);

        // // Perform element-wise division safely
        // Eigen::VectorXd phi_ratio = qrSolverPtr_->phi.array() / phi_solution_segment.array();

        // // Output the result with 10 elements per line
        // std::cout << "from URDF parameter / solution (element-wise):" << std::endl;
        // for (int i = 0; i < phi_ratio.size(); ++i) {
        //     std::cout << phi_ratio(i) << " ";
            
        //     if ((i + 1) % 10 == 0) {
        //         std::cout << std::endl;
        //     }
        // }



    }
    catch (std::exception& e) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    return 0;
}
