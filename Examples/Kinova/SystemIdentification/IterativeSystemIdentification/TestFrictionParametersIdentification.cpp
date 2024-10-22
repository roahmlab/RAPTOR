#include "FrictionParametersIdentification.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

using namespace RAPTOR;

int main(int argc, char* argv[]) {
    // Initialize model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);

    model.friction.setZero();
    model.damping.setZero();
    model.armature.setZero(); 

    // Initialize Regroup matrix
    auto qrSolverPtr_ = std::make_shared<QRDecompositionSolver>(model); 
    qrSolverPtr_->generateRandomObservation(40000);
    qrSolverPtr_->computeRegroupMatrix();
    // const std::string solFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/full_parameters_solution_"
    //                             + std::string(argv[1]) +".csv";
    //     Eigen::VectorXd phi = Utils::initializeEigenMatrixFromFile(solFile).col(0);
    // std::shared_ptr<Eigen::VectorXd> phiPtr_ = std::make_shared<Eigen::VectorXd>(phi);

    // Initialize data
    bool include_offset_input = true ;
    
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

    torqueData = -torqueData;
    std::cout << torqueData.row(0)<<std::endl;


    int N = posData.rows();

    std::shared_ptr<Eigen::MatrixXd> posDataPtr_ = std::make_shared<Eigen::MatrixXd>(posData.transpose());
    std::shared_ptr<Eigen::MatrixXd> velDataPtr_ = std::make_shared<Eigen::MatrixXd>(velData.transpose());
    std::shared_ptr<Eigen::MatrixXd> accDataPtr_ = std::make_shared<Eigen::MatrixXd>(accData.transpose());
    std::shared_ptr<Eigen::MatrixXd> torqueDataPtr_ = std::make_shared<Eigen::MatrixXd>(torqueData.transpose());
    // compute nominal torque from data
    // for (Index i = 0; i < N; i++) {
    //     const Eigen::VectorXd& q = posDataPtr_->col(i);
    //     const Eigen::VectorXd& v = velDataPtr_->col(i);
    //     const Eigen::VectorXd& a = accDataPtr_->col(i);
    //     pinocchio::rnea(model, data, q, v, a);
    //     torqueDataPtr_->col(i) = data.tau + Eigen::VectorXd::Random(model.nv);
    // } 

    // Eigen::MatrixXd tau_inertials = Eigen::MatrixXd::Zero(model.nv, N);
    // for (int i = 0; i < N; i++) {
    //     Eigen::VectorXd q(7);
    //     q << 1.0011758805, 0.0913270190, -1.6488381624, 2.3811290264, 1.8223775625, 0.1466418207, 0.9315121174;
    //     Eigen::VectorXd v(7);
    //     v << 0.0000005651, -0.0000002095, -0.0000001233, -0.0000000360, -0.0000001823, 0.0000000111, 0.0000000027;
    //     Eigen::VectorXd a(7);
    //     a << 0.0018199695, -0.0006747015, -0.0003970416, -0.0001159423, -0.0005872960, 0.0000357368, 0.0000085602;


    //     pinocchio::rnea(model, data, q, v, a);

    //     tau_inertials.col(i) = data.tau;
    //     if (i ==0 or i==1) {
    //         std::cout << (data.tau).transpose() << std::endl;
    //         std::cout << (torqueDataPtr_->col(i)).transpose() << std::endl;

    //     }
    // } 

    // Initialize Kinova optimizer
    SmartPtr<FrictionParametersIdentification> mynlp = new FrictionParametersIdentification();
    try {
	    mynlp->set_parameters(model,
                              posDataPtr_,
                              velDataPtr_,
                              accDataPtr_,
                              torqueDataPtr_,
                            //   qrSolverPtr_,
                              include_offset_input);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
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
    app->Options()->SetStringValue("ma57_automatic_scaling", "yes");
    if (mynlp->enable_hessian) {
        app->Options()->SetStringValue("hessian_approximation", "exact");
    }
    else {
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    }

    // // For gradient checking
    app->Options()->SetStringValue("output_file", "ipopt.out");
    app->Options()->SetStringValue("derivative_test", "second-order");
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

        std::cout << "paramer solution: " << (mynlp->solution).transpose() << std::endl;

        const std::string outputfolder = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/";
        std::ofstream solution(outputfolder + "friction_parameters_solution_" + std::string(argv[1]) + ".csv");

        solution << std::setprecision(16);
        for (int i = 0; i < mynlp->solution.size(); i++) {
            solution << mynlp->solution(i) << std::endl;
        }

        Eigen::VectorXd friction = mynlp->solution.head(model.nv);
        Eigen::VectorXd damping = mynlp->solution.segment(model.nv, model.nv);
        Eigen::VectorXd armature = mynlp->solution.segment(2 * model.nv, model.nv);
        Eigen::VectorXd offset = Eigen::VectorXd::Zero(model.nv);
        if (include_offset_input) {
            offset = mynlp->solution.tail(model.nv);
        }

        // Eigen::MatrixXd FullObservationMatrix = Eigen::MatrixXd::Zero(N * model.nv, 10 * model.nv);
        // for (int i = 0; i < N; i++) {
        //     const Eigen::VectorXd& q_d = velDataPtr_->col(i);
        //     const Eigen::VectorXd& q_dd = accDataPtr_->col(i);
        //     const Eigen::VectorXd& q = posDataPtr_->col(i);

        //     pinocchio::computeJointTorqueRegressor(
        //         model, data,
        //         q, q_d, q_dd);

        //     FullObservationMatrix.middleRows(i * model.nv, model.nv) = 
        //         data.jointTorqueRegressor;
        // }

        // Eigen::MatrixXd RegroupedObservationMatrix = 
        // FullObservationMatrix * qrSolverPtr_->RegroupMatrix;

        // Eigen::VectorXd tau_inertials = Eigen::VectorXd::Zero(10 * model.nv);
    
        // tau_inertials = RegroupedObservationMatrix * qrSolverPtr_->beta;
   


        // Eigen::VectorXd phi = Eigen::VectorXd::Zero(10 * model.nv);
        // for (int i = 0; i < model.nv; i++) {
        //     const int pinocchio_joint_id = i + 1;
        //     phi.segment<10>(10 * i) = 
        //         model.inertias[pinocchio_joint_id]
        //             .toDynamicParameters();
        // }
        // Eigen::VectorXd tau_inertials = Eigen::VectorXd::Zero(10 * model.nv);
    
        // tau_inertials = FullObservationMatrix * phi;
   

        Eigen::MatrixXd tau_inertials = Eigen::MatrixXd::Zero(model.nv, N);
        for (int i = 0; i < N; i++) {
            const Eigen::VectorXd& q = posDataPtr_ ->col(i);
            const Eigen::VectorXd& v = velDataPtr_->col(i);
            const Eigen::VectorXd& a = accDataPtr_->col(i);

            pinocchio::rnea(model, data, q, v, a);

            tau_inertials.col(i) = data.tau;
        } 


        const std::string outputfolder1 = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/";
        std::ofstream estimate_tau(outputfolder1 + "friction_estimate_tau_"  + std::string(argv[1]) +".csv");

        for (Index i =0 ; i < N; i++) {
            const Eigen::VectorXd& q_d = velDataPtr_->col(i);
            const Eigen::VectorXd& q_dd = accDataPtr_->col(i);
            const Eigen::VectorXd& tau = torqueDataPtr_->col(i);
            const Eigen::VectorXd& q = posDataPtr_->col(i);
            const Eigen::VectorXd& tau_inertial = tau_inertials.col(i);
            // const Eigen::VectorXd& tau_inertial = tau_inertials.segment(i * model.nv, model.nv);
            Eigen::VectorXd tau_estimated = Eigen::VectorXd::Zero(model.nv);

            Eigen::VectorXd total_friction_force = 
                (friction.cwiseProduct(q_d.cwiseSign()) +
                damping.cwiseProduct(q_d) +
                armature.cwiseProduct(q_dd) +
                offset);
            tau_estimated = tau_inertial + total_friction_force;

            for (int j = 0; j < tau_estimated.size(); ++j) {
                estimate_tau << tau_estimated(j);
                if (j < tau_estimated.size() - 1) {
                    estimate_tau << " ";  
                }
            }
            estimate_tau << std::endl;
        }

    // Eigen::VectorXd tau_inertials1 = Eigen::VectorXd::Zero(model.nv);

    // for (int i = 0; i < N; i++) {
    //     Eigen::VectorXd q(7);
    //     q << 1.00118076, 0.09136623,-1.64880621, 2.38126,1.82228,0.14661227,0.93155;
    //     // q= Eigen::VectorXd::Zero(7);
    //     Eigen::VectorXd v(7);
    //     // v << 0.403876417605741, -0.183062529759288, -0.0943837857317314, -8.15005777823162e-05, -0.124105515843069, 0.0102116544925587, -0.00524206000582545;
    //     v = Eigen::VectorXd::Zero(7);
    //     Eigen::VectorXd a(7);
    //     // a << -1.92276862278343, 1.52821104080408, -3.33140487520136, 0.11742782389603, -2.8557732836922, 0.253655612030757, 0.126359311385054;
    //     a = Eigen::VectorXd::Zero(7);
      
    //     pinocchio::rnea(model, data, q, v, a);

    //     tau_inertials1 = data.tau;

    //     if (i ==0 ) {
    //         std::cout << tau_inertials1.transpose() << std::endl;

    //     }

    //     // Eigen::VectorXd total_friction_force = 
    //     //         (friction.cwiseProduct(v.cwiseSign()) +
    //     //         damping.cwiseProduct(v) +
    //     //         armature.cwiseProduct(a) +
    //     //         offset);

    //     // Eigen::VectorXd tau_estimated = tau_inertials1 + total_friction_force;
    //     // if (i ==0 ) {
    //     //     std::cout << tau_estimated.transpose() << std::endl;

    //     // }
    // } 
    
        
    }
    catch (std::exception& e) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    return 0;
}
