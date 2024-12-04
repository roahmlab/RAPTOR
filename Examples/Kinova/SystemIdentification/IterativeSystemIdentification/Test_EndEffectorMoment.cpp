#include "MomentumRegressor.h"
#include "TrajectoryData.h"
#include "Plain.h"
#include "pinocchio/algorithm/rnea.hpp"
#include "EndEffectorMoment.h"
#include <omp.h>

using namespace RAPTOR;

// Function declarations
Eigen::MatrixXd theta_to_LMI(Eigen::VectorXd& theta);
double BuresWassersteinDistance(Eigen::MatrixXd& A, Eigen::MatrixXd& B);
Eigen::MatrixXd matrixSquareRoot(Eigen::MatrixXd& A);
Eigen::VectorXd x_to_theta(Eigen::VectorXd& x);

int main(int argc, char* argv[]) {
    // check if the file number is provided
    if (argc < 2) {
            std::cerr << "Error: No arguments provided. please choose the downsampled file" << std::endl;
            return 1; // 
        }

    std::string file_number = std::string(argv[1]);


    // Load the robot model
    const std::string urdf_filename  = "/workspaces/RAPTOR/Robots/kinova-gen3/kinova_grasp_fixed.urdf";

    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);
    pinocchio::Data data(model);

    // Disable rotor inertia, friction, and damping
    model.friction.setZero();
    model.damping.setZero();
    model.armature.setZero();

    Eigen::VectorXd phi = Eigen::VectorXd::Zero(10 * model.nv);
    for (int i = 0; i < model.nv; i++) {
        const int pinocchio_joint_id = i + 1;
        phi.segment<10>(10 * i) =
            model.inertias[pinocchio_joint_id]
                .toDynamicParameters();
    }


    std::ofstream inertia_parameters("inertia_parameters.csv");
    for (int i = 0; i < phi.size(); i++) {
        inertia_parameters << phi(i) << std::endl;
    }

    // load the data
    std::string filename = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_effector_params_data/2024_11_20_gripper_id_";
    std::shared_ptr<TrajectoryData> trajPtr_ =
        std::make_shared<TrajectoryData>(filename + file_number + "_together.txt");

    Eigen::MatrixXd ctrl = Utils::initializeEigenMatrixFromFile(filename + file_number + "_together.txt").rightCols(model.nv);

    // Z is not used in the compute function chosen randomly
    std::srand(std::time(nullptr));
    Eigen::VectorXd z = M_2_PI * Eigen::VectorXd::Random(7);  // meaningless
    
    // load  friction parameters
    const std::string solFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_moment_id_data/friction_parms.csv";
    Eigen::VectorXd fricition_parameters = Utils::initializeEigenMatrixFromFile(solFile).col(0);
    std::cout << "friction_parameters size: " << fricition_parameters.size() << std::endl;

    if (fricition_parameters.size() > 4 * model.nv) {
            std::cerr << " Friction solution file is wrong" << std::endl;
            return 1; // 
    }

    bool include_offset_input = (fricition_parameters.size() == 4 * model.nv)? true : false;
    std::cout << "include_offset_input: " << include_offset_input << std::endl; 

    const Eigen::VectorXd friction = fricition_parameters.head(model.nv);
    const Eigen::VectorXd damping = fricition_parameters.segment(model.nv, model.nv);
    const Eigen::VectorXd armature  = fricition_parameters.segment(2 * model.nv, model.nv);
    Eigen::VectorXd offset = Eigen::VectorXd::Zero(model.nv);
    if (include_offset_input) {
        offset  = fricition_parameters.tail(model.nv);
    }

    // Initialize MomentumRegressor
    std::shared_ptr<MomentumRegressor> ridPtr_ = std::make_shared<MomentumRegressor>(model, trajPtr_);

    // Compute momentums regressor
    ridPtr_ -> compute(z, false);

    // gravity regressor
    std::shared_ptr<TrajectoryData> trajPtr2_ =
        std::make_shared<TrajectoryData>(filename + file_number + "_together.txt");
    for (int i = 0; i < trajPtr2_->N; i++) {
        trajPtr2_->q_d(i).setZero();
    }
    std::shared_ptr<RegressorInverseDynamics> Y_gPtr_ = std::make_shared<RegressorInverseDynamics>(model, 
                                                            trajPtr2_,
                                                            false); 
    Y_gPtr_ -> compute(z, false);
   
    // Choose the integration horizion
    int H = 5;

    // Integration horizon
    int num_seg = (trajPtr_->N) / H;

    if (num_seg <= 0) {
        std::cerr << "Error: num_seg is zero or negative. Check H and trajPtr_->N." << std::endl;
        return 1;
    }

    // Initialize matrices and vectors
    Eigen::MatrixXd A_full = Eigen::MatrixXd::Zero(model.nv * num_seg, 10 * model.nv);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(model.nv * num_seg);
   
    // Create cell-like structures using std::vector
    std::vector<Eigen::MatrixXd> Aelt(num_seg);
    std::vector<Eigen::VectorXd> Belt(num_seg);

    int i = 0;
    #pragma omp parallel for shared(trajPtr_,  ridPtr_ , A_full , Y_gPtr_ ,b, Aelt, Belt, ctrl, friction, damping, armature) private(i) schedule(dynamic, 1)
    for (i = 0; i < trajPtr_->N - H-1; i += H) {
        int seg_start = i;
        int seg_end = seg_start + H;

        Eigen::MatrixXd Y_Hqd_1 = ridPtr_ -> Y.middleRows(seg_start * model.nv, model.nv);
        Eigen::MatrixXd Y_Hqd_2 = ridPtr_ -> Y.middleRows(seg_end * model.nv, model.nv);

        Eigen::MatrixXd int_Y_CTqd_g = Eigen::MatrixXd::Zero(model.nv, 10 * model.nv);
        Eigen::VectorXd int_ctrl = Eigen::VectorXd::Zero(model.nv);


        for (int j = seg_start; j< seg_end; j++) {
            double dt = trajPtr_->tspan(j+1) - trajPtr_->tspan(j);

            Eigen::MatrixXd Y_CTqd_i = ridPtr_ -> Y_CTv.middleRows(j * model.nv, model.nv);
            Eigen::MatrixXd Yg_i = Y_gPtr_ ->Y.middleRows(j * model.nv, model.nv);

            int_Y_CTqd_g += (Y_CTqd_i - Yg_i) * dt;

            int_ctrl += ctrl.row(j).transpose() * dt
                        - friction.cwiseProduct(trajPtr_->q_d(j).cwiseSign()) * dt
                        - damping.cwiseProduct(trajPtr_->q_d(j)) * dt;

            if (include_offset_input) {
                int_ctrl -= offset * dt;
            }
        }
        int s = i / H;
        Aelt[s] = (Y_Hqd_2 - Y_Hqd_1) - int_Y_CTqd_g;
        Belt[s] = int_ctrl - armature.cwiseProduct(trajPtr_->q_d(seg_end) - trajPtr_->q_d(seg_start));
    }

    // Assemble the global matrices A_full, A, and b
    for (int s = 0; s < num_seg-1; ++s) { // check later, the aelt(s) at end is 
        A_full.middleRows(s * model.nv, model.nv) = Aelt[s];
        b.segment(s * model.nv, model.nv) = Belt[s];
    }

    // Write the matrices to file
    std::string filename_out = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_moment_id_data/gripper_" + file_number;
    Utils::writeEigenMatrixToFile(A_full.middleRows(0, (num_seg-1)* model.nv), filename_out + "_A_full_H_5.csv");
    Eigen::MatrixXd b_matrix = Eigen::Map<Eigen::MatrixXd>(b.data(), (num_seg-1)* model.nv, 1);
    Utils::writeEigenMatrixToFile(b_matrix, filename_out + "_b_H_5.csv");

    // laod the file, and can comment the above code
    // bool include_offset_input = false;
    // Eigen::MatrixXd A_full = RAPTOR::Utils::initializeEigenMatrixFromFile("A_full.csv");
    // Eigen::VectorXd b = RAPTOR::Utils::initializeEigenMatrixFromFile("b.csv").leftCols(1);
    // Eigen::VectorXd phi = RAPTOR::Utils::initializeEigenMatrixFromFile("inertia_parameters.csv").leftCols(1);

    std::shared_ptr<Eigen::VectorXd> inertia_parametersPtr_ = std::make_shared<Eigen::VectorXd>(phi);
    std::shared_ptr<Eigen::VectorXd> b_ = std::make_shared<Eigen::VectorXd>(b);
    std::shared_ptr<Eigen::MatrixXd> A_full_ = std::make_shared<Eigen::MatrixXd>(A_full);
  
    // Initialize the Ipopt problem
    SmartPtr<EndEffectorMoment> mynlp = new EndEffectorMoment();
    try {
	    mynlp->set_parameters(A_full_,
                              b_,
                              inertia_parametersPtr_,
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
	app->Options()->SetNumericValue("max_wall_time", 100);
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
    // // app->Options()->SetStringValue("derivative_test", "second-order");
    // app->Options()->SetStringValue("derivative_test", "first-order");
    // app->Options()->SetNumericValue("derivative_test_perturbation", 1e-8);
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-5);
    // app->Options()->SetNumericValue("point_perturbation_radius", 1);

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
        std::cout <<"solution x" << mynlp->solution.transpose()<< std::endl;
        std::cout << "paramer solution: " << x_to_theta(mynlp->solution).transpose() << std::endl;
       


        Eigen::VectorXd phi_tail = phi.tail(10);
        Eigen::MatrixXd LMI_original = theta_to_LMI(phi_tail);
        Eigen::VectorXd theta_sol = x_to_theta(mynlp->solution);
        Eigen::MatrixXd LMI_solution = theta_to_LMI(theta_sol);
        double distance = BuresWassersteinDistance(LMI_original, LMI_solution);
        std::cout << "Bures-Wasserstein distance: " << distance << std::endl;

        phi.tail(10) = theta_sol;
        // std::cout << "residual" << (A_full * phi - b) << std::endl;                              
        std::cout << "average residual" << (A_full * phi - b).cwiseAbs().sum()/b.size()<< std::endl;
        std::cout << "max residual" << (A_full * phi - b).cwiseAbs().maxCoeff() << std::endl;

        
    }
    catch (std::exception& e) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    return 0;

}

// helfer functions
double BuresWassersteinDistance(Eigen::MatrixXd& A, Eigen::MatrixXd& B) {
    if (A.rows() != A.cols() || B.rows() != B.cols() || A.rows() != B.rows()) {
        throw std::invalid_argument("Matrices A and B must be square and of the same size.");
    }

    Eigen::MatrixXd sqrtA = matrixSquareRoot(A);
    Eigen::MatrixXd middle = sqrtA * B * sqrtA;
    Eigen::MatrixXd sqrtMiddle = matrixSquareRoot(middle);

    double traceA = A.trace();
    double traceB = B.trace();
    double traceSqrtMiddle = sqrtMiddle.trace();

    return std::sqrt(traceA + traceB - 2 * traceSqrtMiddle);
}

Eigen::MatrixXd matrixSquareRoot(Eigen::MatrixXd& matrix) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();

    // Compute the square root of the diagonal matrix of eigenvalues
    Eigen::VectorXd sqrtEigenvalues = eigenvalues.array().sqrt();

    // Reconstruct the square root matrix
    return eigenvectors * sqrtEigenvalues.asDiagonal() * eigenvectors.transpose();
}

Eigen::MatrixXd theta_to_LMI(Eigen::VectorXd& theta) {
    // Extract mass
    double mass = theta(0);

    // Extract center of mass (com)
    Eigen::Vector3d com = theta.segment<3>(1);

    // Extract inertia elements
    double Ixx = theta(4);
    double Ixy = theta(5);
    double Iyy = theta(6);
    double Ixz = theta(7);
    double Iyz = theta(8);
    double Izz = theta(9);

    // Construct the inertia matrix
    Eigen::Matrix3d inertia;
    inertia << Ixx, Ixy, Ixz,
               Ixy, Iyy, Iyz,
               Ixz, Iyz, Izz;

    // Compute trace of the inertia matrix
    double inertia_trace = inertia.trace();

    // Compute the LMI matrix
    Eigen::MatrixXd LMI(4, 4);
    LMI.setZero();
    LMI.block<3, 3>(0, 0) = 0.5 * inertia_trace * Eigen::Matrix3d::Identity() - inertia;
    LMI.block<3, 1>(0, 3) = com;                                                        
    LMI.block<1, 3>(3, 0) = com.transpose();                                          
    LMI(3, 3) = mass;                                                                  

    return LMI;
}

Eigen::VectorXd x_to_theta(Eigen::VectorXd& z) {

    double d1 = z[0];
    double d2 = z[1];
    double d3 = z[2];
    double d4 = z[3];
    double s12 = z[4];
    double s23 = z[5];
    double s13 = z[6];
    double t1 = z[7];
    double t2 = z[8];
    double t3 = z[9];

    Eigen::Matrix4d U;
    U << std::exp(d1), s12,      s13,      t1,
        0.0,     std::exp(d2),  s23,      t2,
        0.0,     0.0,      std::exp(d3),  t3,
        0.0,     0.0,      0.0,      std::exp(d4);

    // Compute LMI = U' * U
    Eigen::Matrix4d LMI = U.transpose() * U;

    // End-effector parameters
    Eigen::VectorXd theta = Eigen::VectorXd::Zero(10);
    theta(0) = LMI(3, 3);                        
    theta.segment<3>(1) = LMI.block<3, 1>(0, 3); 
    theta(4) = LMI(1, 1) + LMI(2, 2);            // IXX
    theta(5) = -LMI(0, 1);                       // IXY
    theta(6) = LMI(0, 0) + LMI(2, 2);            // IYY
    theta(7) = -LMI(0, 2);                       // IXZ
    theta(8) = -LMI(1, 2);                       // IYZ
    theta(9) = LMI(0, 0) + LMI(1, 1);            // IZZ

    return theta;
}
