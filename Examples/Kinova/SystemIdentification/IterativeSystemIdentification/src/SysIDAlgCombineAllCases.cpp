#include <vector>
#include <Eigen/Dense>
#include "Utils.h"
#include <cmath>
#include <iostream> 
#include <memory>
#include <cstdio>
#include <cstdlib>

#include "BaseParametersIdentification.h"
#include "FrictionParametersIdentification.h"
#include "SysIDAlgCombineAllCases.h"


namespace RAPTOR {
namespace Kinova {
/* ---------------------------------- input --------------------------------- */
/*
N: number of data points
Alg_case: 0 is full links, 1 is end_effector, 2 is friction only
error_threshold: hard redescender threshold
Weight_tol: weight tolerate, stop criterion in the first while loop of baseparametersidentification 
O_sqrt_tol: covariance tolerate, stop criterion in the second while loop of baseparametersidentification
Alpha_tol: stop criterion in the  frictionparametersidentification
include_friction_offset: wether the friction has offet or not 
files: [p, v, a , solution]  // later will add Torque from experiment 
*/

SysIDAlgCombineAllCases::SysIDAlgCombineAllCases(
        const int N,
        const int Alg_case, 
        const double error_threshold, //hard redescender
        const double Weight_tol,
        const double Alpha_tol,
        const double O_sqrt_tol,
        const bool include_friction_offset,
        std::vector<std::string> files,
        Model &model_input
)
    :   N_(N), Alg_case_(Alg_case), 
        error_threshold_(error_threshold), Weight_tol_(Weight_tol), Alpha_tol_(Alpha_tol),O_sqrt_tol_(O_sqrt_tol),
        include_friction_offset_(include_friction_offset)
{
    
    /* --------------  load the optimized trajectory information -------------- */
    position= files[0];
    velocity = files[1];
    acceleration = files[2];
    solution = files[3];

    /* -------------------------  regroup informaiton ------------------------- */
    /*
    pi_d: the vector of depentent parameters
    Ginv: inverse of bijective map to get inertial parameters from base inertial parameters
    Ad: selection matrix for dependant and unidentifiable parameters
    Aid: selection matrix for indepenent parameters
    Kd: regrouping transformation matrix
    b_dim: the number of independent parameters
    d_dim: the number of dependent parameters
    */
    regroupPtr_ = std::make_shared<QRDecompositionSolver>();
    if (Alg_case_ ==0){
        regroupPtr_ ->getData();
        pi_d_ = regroupPtr_->pi_d;
        Ginv_ = regroupPtr_->Ginv;
        Aid_ = regroupPtr_->Aid;
        Ad_ = regroupPtr_->Ad;
        Kd_ = regroupPtr_->Kd;
        b_dim_ = regroupPtr_->dim_id;
        d_dim_ = regroupPtr_->dim_d;
    }

    /* ----------------------------- set parameters ----------------------------- */
    // set parameters for full inertia, end_effector and friction case
    /*
    nLinks: number of joint in the robot 
    p_ip: number of inertia parameters
    b_full: number of regroup inertia parameters + number of motor parameters
    p_full: number of inertia parameters + number of motor parameters
    Wfull: observation matrix for inertial, motor, friction and offset parameters
    alphanew: new friction exponents from optimization
    */
    fm_num_ = include_friction_offset_? 4:3;
    modelPtr_ = std::make_shared<Model>(model_input);
    switch (Alg_case_) {
        case 0:
            nLinks_= modelPtr_->nv;
            p_ip_ =10 * nLinks_;
            b_full_ = b_dim_ + fm_num_ * nLinks_;
            p_full_ = p_ip_ + fm_num_ * nLinks_;
            Wfull_ = MatX::Zero(modelPtr_->nv *N_, b_full_);
            alphanew_ = VecX::Ones(nLinks_);
            break;
        case 1:
            nLinks_= 1;
            p_ip_ = 10;
            b_full_ = p_ip_+ fm_num_ * nLinks_;
            p_full_ = p_ip_ + fm_num_ * nLinks_;
            Wfull_ = MatX::Zero(modelPtr_->nv *N_, b_full_);
            alphanew_ = VecX::Ones(nLinks_);
            break;
        case 2:
            nLinks_= 1;
            p_ip_ = 0;
            b_full_ = p_ip_+ fm_num_ * nLinks_;
            p_full_ = p_ip_ + fm_num_ * nLinks_;
            Wfull_ = MatX::Zero(modelPtr_->nv *N_, b_full_);
            alphanew_ = VecX::Ones(nLinks_);
            break;
        default:
            throw std::runtime_error("Invalid mode.");
    }

  
    /* ------------  Initialize the full parameters from the model ------------ */

    // Structure [Inertia*nLinks_,Ia*nLinks_, firction *nLinks_]
    X0_1_ = VecX::Zero(p_full_);

    // intertia
    if (Alg_case_ == 0){
        for (int i = 0; i < nLinks_; i++) {
            X0_1_.segment<10>(10 * i) = model_input.inertias[i + 1].toDynamicParameters();
        }
    }else if (Alg_case_ == 1) {
        X0_1_.segment<10>(0) = model_input.inertias[nLinks_].toDynamicParameters();
    }

    // motor parameters
    for (int i = 0; i < nLinks_; i++) {
        // Ia
        X0_1_.segment<1>(p_ip_+ i) << modelPtr_->rotorInertia(i);
        // friction 
        if (include_friction_offset_){
            X0_1_.segment(p_ip_ + nLinks_ + (fm_num_ - 1) * i, fm_num_ -1) << modelPtr_->friction(i), modelPtr_->damping(i), 0.0;
        }else {
            X0_1_.segment(p_ip_ + nLinks_ + (fm_num_ -1) * i, fm_num_ -1) << modelPtr_->friction(i), modelPtr_->damping(i);
        }
    }


    /* ------------------------- lower and upper bounds ------------------------- */
    
    lb_ = VecX::Zero(b_full_);
    ub_ = VecX::Zero(b_full_);
 
    if (Alg_case_ == 0){
        // Initialize lower and upper bounds, after regroup
        VecX lowerB;
        VecX upperB;
        // Inertial parameters ±10% of initial value
        lowerB =X0_1_.segment(0,p_ip_) * 0.9;
        upperB = X0_1_.segment(0,p_ip_) * 1.1;

        // Regroup the bounds for inertial parameters
        lb_.segment(0, b_dim_) = Aid_.transpose() * lowerB + Kd_ * Ad_.transpose() * lowerB;
        ub_.segment(0, b_dim_) = Aid_.transpose() * upperB + Kd_ * Ad_.transpose() * upperB;
    

        // Check and switch regrouped bounds 
        for (int j = 0; j < b_dim_; ++j) {
            if (lb_(j) >= ub_(j)) {
                std::swap(lb_(j), ub_(j));
            }
        }
        // motor parameters bound [0,20]
        lb_.segment(b_dim_, fm_num_ * nLinks_).setConstant(0.0);
        ub_.segment(b_dim_, fm_num_ * nLinks_).setConstant(20.0);

    }else{
        // Inertial parameters ±10% of initial value
        lb_.segment(0, p_ip_) =X0_1_.segment(0, p_ip_) * 0.9;
        ub_.segment(0, p_ip_)  = X0_1_.segment(0, p_ip_) * 1.1;

        for (int j = 0; j < p_ip_; ++j) {
            if (lb_(j) >= ub_(j)) {
                std::swap(lb_(j), ub_(j));
            }
        }

        // motor parameters bound [0,20]
        lb_.segment( p_ip_, fm_num_ * nLinks_).setConstant(0.0);
        ub_.segment( p_ip_, fm_num_ * nLinks_).setConstant(20.0);
    }

    /* ------------------ process Observation matrix and Torque ----------------- */
    ridPtr_ = std::make_shared<RegressorInverseDynamics>(model_input,
                                                        position,
                                                        velocity,
                                                        acceleration);

    // ineratia part of Observation matrix ( only full and end_effector have)
    ridPtr_->computeWithDataImport(solution);

    // motor part of Observation matrix
    MatX Wfm, Wm, Wf; 
    Wextension(alphanew_, Wfm, Wm, Wf);

    switch (Alg_case_) {
        case 0:
            W_ip_ =ridPtr_->Y.leftCols(10*nLinks_);
            // T_ = W_ip_* X0_1_.head(10*nLinks_);
            regroupParameters();
            Wfull_ << Wb_, Wfm;
            T_ =  Wfull_ * X0_1_;   // should read from machine 
            break;
        case 1:
            W_ip_ =ridPtr_->Y.middleCols((modelPtr_->nv - 1) * 10, 10);
            // T_ = W_ip_* X0_1_.head(10);
            Wfull_ << W_ip_ , Wfm;
            T_ =  Wfull_ * X0_1_;  // should read from machine 
            break;
        case 2:
            Wfull_ << Wfm;
            T_ =  Wfull_ * X0_1_;  // should read from machine 
            break;
    }
    
}

void SysIDAlgCombineAllCases::runAlgorithm() {
    //  initialize alpha
    VecX alphaold = VecX::Zero(nLinks_);
    VecX alphanew = alphanew_; 

    while ((alphaold -alphanew).norm() > Alpha_tol_) {
        std::cout<<"(alphaold -alphanew).norm()"<< (alphaold -alphanew).norm()<<std::endl;
        // first identtifcation weight and base parameters
        runFirstOptimization();
        
        // second identifcation frictionexponent
        runSecondOptimization(alphaold, alphanew);
        // break;
    }
}

void SysIDAlgCombineAllCases::regroupParameters() {

    // regroup w_ip
    Wb_ = W_ip_ * Aid_;

    // updata params
    VecX X0_ip = Aid_.transpose() * X0_1_.head(10 * nLinks_) + Kd_ * Ad_.transpose() * X0_1_.head(10 * nLinks_);
    VecX X0_fm = X0_1_.tail(X0_1_.size() - 10 * nLinks_);

    X0_1_.resize(b_full_);
    X0_1_.head(b_dim_) = X0_ip;
    X0_1_.tail(X0_fm.size()) = X0_fm;
}

void SysIDAlgCombineAllCases::Wextension(const VecX& alpha, MatX& Wfm, MatX& Wm, MatX& Wf) {
    int wf_cols = include_friction_offset_? 3 * nLinks_ :2 * nLinks_;
    Wf = MatX::Zero(modelPtr_->nv * N_, wf_cols);
    Wm = MatX::Zero(modelPtr_->nv * N_, nLinks_);

    for (int i = 0; i < N_; ++i) {
        const VecX& q = ridPtr_->trajPtr_q(i);
        const VecX& q_d = ridPtr_->trajPtr_q_d(i);
        const VecX& q_dd = ridPtr_->trajPtr_q_dd(i);
        for (int j =0; j< nLinks_ ; ++j){
            Wm(i * modelPtr_->nv + j, j) = q_dd(j); 
            if (include_friction_offset_) {
                Wf(i * modelPtr_->nv + j, 3 * j) = Utils::sign(q_d(j));
                Wf(i * modelPtr_->nv + j, 3 * j + 1) = std::pow(std::abs(q_d(j)), alpha(j));
                Wf(i * modelPtr_->nv + j, 3 * j + 2) = 1;
            }else{
                Wf(i * modelPtr_->nv + j, 3 * j) = Utils::sign(q_d(j));
                Wf(i * modelPtr_->nv + j, 3 * j + 1) = std::pow(std::abs(q_d(j)), alpha(j));
            }
        }
    }
    Wfm.resize(Wm.rows(), Wm.cols() + Wf.cols());
    Wfm << Wm, Wf;
}

/* -------------------------- runFirstOptimization -------------------------- */
/*
Dnew, Dold: the weight matrix for the toruqe
De : the wieght matrix for the regressor matrix
Oold, Onew : error vector covarianace matrix
Rh: normialized residual 
Eh: reshaped normialized residual (nLinks_ * N_)
Ts: normialized toruqe
Ws: normialized regressor matrix
*/
void SysIDAlgCombineAllCases::runFirstOptimization() {
    // Initialize
    VecX Dnew;
    VecX Dold;
    VecX Rh;
    MatX Eh;
    MatX De;
    MatX Oold ;
    MatX Onew ;
    VecX Ts;
    MatX Ws;

    Dnew = VecX::Ones(N_*nLinks_);
    Dold = VecX::Zero(N_*nLinks_);

    De = Dnew.replicate(1,b_full_);

    // start the loop for the weight
    while ((Dold - Dnew).norm() > Weight_tol_){
        std::cout<<"(Dold - Dnew).norm() "<< (Dold - Dnew).norm()<<std::endl;
        //  Initialize covarianace
        Oold = MatX::Zero(nLinks_ , nLinks_ );
        Onew = MatX::Identity(nLinks_ , nLinks_ );  
        MatX O_sqrt = Onew.llt().matrixL();

        // start loop for covariance matrix
        while ((Onew - Oold).norm() > O_sqrt_tol_){
            std::cout<<"(Onew - Oold).norm()"<< (Onew - Oold).norm()<<std::endl;
            // normialized toruqe and regressor matrix
            Ts = VecX::Zero(nLinks_*N_);
            Ws = MatX::Zero(nLinks_ * N_, Wfull_.cols());  //move above
            for (int i =0; i < N_; ++i){
                Ts.segment(i*nLinks_, nLinks_) = O_sqrt.inverse()* T_.segment(i*nLinks_,nLinks_);
                Ws.block(i*nLinks_, 0, nLinks_, Wfull_.cols()) = O_sqrt.inverse() * Wfull_.block(i * nLinks_, 0, nLinks_, Wfull_.cols());
            }
        
            //  torque and regressor after weight 
            Th_ = Dnew.array() * Ts.array();
            Wh_ = De.array() * Ws.array(); 

            //optimise base parameters  
            SmartPtr<BaseParametersIdentification> mynlp = new BaseParametersIdentification();
            try {
                mynlp->BaseParametersIdentification::set_parameters(
                                    Wh_,
                                    Th_,
                                    X0_1_,
                                    include_friction_offset_,
                                    *modelPtr_,
                                    regroupPtr_,
                                    lb_,
                                    ub_,
                                    b_full_,
                                    fm_num_,
                                    Alg_case_);
                mynlp->constr_viol_tol = 1e-5;
            }
            catch (std::exception& e) {
                std::cerr << e.what() << std::endl;
                throw std::runtime_error("BaseParametersIdentification: Error initializing Ipopt class! Check previous error message!");
            }

            SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
            app->Options()->SetNumericValue("tol", 1e-6);
            app->Options()->SetNumericValue("constr_viol_tol", mynlp->constr_viol_tol);
            app->Options()->SetNumericValue("max_wall_time", 60);
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
            app->Options()->SetStringValue("output_file", "ipopt.out");
            app->Options()->SetStringValue("derivative_test", "first-order");
            app->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
            app->Options()->SetNumericValue("derivative_test_tol", 1e-3);
            app->Options()->SetNumericValue("point_perturbation_radius", 0);

            // Initialize the IpoptApplication and process the options
            ApplicationReturnStatus status;
            status = app->Initialize();
            if( status != Solve_Succeeded ) {
                throw std::runtime_error("Error during initialization of optimization!");
            }

            try {
                status = app->OptimizeTNLP(mynlp);
            }
            catch (std::exception& e) {
                throw std::runtime_error("Error solving optimization problem! Check previous error message!");
            }

            // updata base parameters
            VecX X_sol = mynlp->solution;
            X0_1_ = X_sol;
      
            // calculate normlized residual and reshape into (nlinks * N_)
            Rh= Th_ - Wh_* X0_1_;
            Eh = Eigen::Map<MatX>(Rh.data(), nLinks_, N_);
        
            // update covariance matrix
            Oold = Onew;
            Onew = (Oold *Eh *Eh.transpose() * Oold)/(N_ - b_full_);
          
            // update sqrt covariance matrix
            O_sqrt=  Onew.llt().matrixL();

            if (!O_sqrt.allFinite()) {
                break;
            }
            // break;
        }
        // update weight 
        Dold = Dnew;
        // udpata Dnew by min(0, T clss hard redescender)
        VecX Psi = (Rh.array().abs()<=error_threshold_).select(VecX::Ones(Rh.size()), VecX::Zero(Rh.size()));
        Dnew = Dold.array().min(Psi.array());
        // updata De
        De = Dnew.replicate(1, b_full_);
        // break;
    }
}

void SysIDAlgCombineAllCases::runSecondOptimization(VecX& alphaold, VecX& alphanew) {
    // X0_2 is X0_1_ with friction parameters zero
    VecX X0_2= VecX::Zero(b_full_);
    int segmentSize = Wfull_.cols()-(fm_num_-1)*nLinks_;
    X0_2.head(segmentSize) = X0_1_.head(segmentSize);

    VecX Fest = Th_ - Wh_ * X0_2;

    //  Xf structure [alpha, friction parameters]
    VecX Xf = X0_1_.tail(fm_num_ * nLinks_);
    Xf.head(nLinks_) = alphanew;

    //optimise funciton 
    SmartPtr<FrictionParametersIdentification> mynlp_second = new FrictionParametersIdentification();
    try {
        mynlp_second->FrictionParametersIdentification::set_parameters(
                                            Xf,
                                            nLinks_,
                                            Fest,
                                            include_friction_offset_,
                                            N_,
                                            ridPtr_ );                   
       mynlp_second->constr_viol_tol = 1e-5;
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("BaseParametersIdentification: Error initializing Ipopt class! Check previous error message!");
    }

    SmartPtr<IpoptApplication> app_second = IpoptApplicationFactory();
    app_second->Options()->SetNumericValue("tol", 1e-6);
    app_second->Options()->SetNumericValue("constr_viol_tol", mynlp_second->constr_viol_tol);
    app_second->Options()->SetNumericValue("max_wall_time", 60.0);
    app_second->Options()->SetIntegerValue("print_level", 5);
    app_second->Options()->SetStringValue("mu_strategy", "adaptive");
    app_second->Options()->SetStringValue("linear_solver", "ma57");
    app_second->Options()->SetStringValue("ma57_automatic_scaling", "yes");
    if (mynlp_second->enable_hessian) {
        app_second->Options()->SetStringValue("hessian_approximation", "exact");
    }
    else {
        app_second->Options()->SetStringValue("hessian_approximation", "limited-memory");
    }

    // For gradient checking
    app_second->Options()->SetStringValue("output_file", "ipopt.out");
    app_second->Options()->SetStringValue("derivative_test", "first-order");
    app_second->Options()->SetNumericValue("derivative_test_perturbation", 1e-7);
    app_second->Options()->SetNumericValue("derivative_test_tol", 1e-3);

    // Initialize the IpoptApplication and process the options
    ApplicationReturnStatus status_second;
    status_second = app_second->Initialize();
    if( status_second != Solve_Succeeded ) {
        throw std::runtime_error("Error during initialization of optimization!");
    }

    try {
        status_second = app_second->OptimizeTNLP(mynlp_second);
    }
    catch (std::exception& e) {
        throw std::runtime_error("Error solving optimization problem! Check previous error message!");
    }

    Xf = mynlp_second->solution;
    // update alpha
    alphaold = alphanew;
    alphanew = Xf.head(nLinks_);
    alphanew_= alphanew;

    // update parameters
    X0_1_.tail((fm_num_-1) * nLinks_) = Xf.tail((fm_num_-1) * nLinks_);
}

Eigen::VectorXd SysIDAlgCombineAllCases::getfinalparam() const {
    return X0_1_;
}

Eigen::VectorXd SysIDAlgCombineAllCases::getAlphaNew() const {
    return alphanew_;
}


} // namespace Kinova
} // namespace RAPTOR
