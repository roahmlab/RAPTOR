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
#include "SysIDAlgFull.h"


namespace RAPTOR {
namespace Kinova {


SysIDAlgFull::SysIDAlgFull(
        // const VecX& lb, //lower bound on paramaters, p_full vector
        // const VecX& ub, //upper bound on paramaters, p_full vector
        // const VecX& lba, //lower bound on friction exponents, n vector
        // const VecX& uba, //upper bound on friction exponents, n vector
        // const VecX& X0_1, //Initial parameter estimation, p_full
        // const VecX& T, //Torque vector. n * m
        const VecX& pi_d,
        // const MatX& W_ip, // Observation matrix of the inertia parameters
        const MatX& Ginv, // The inverse of the bijective mapping that converts base inertia parameters into full inertia parameters
        const MatX& Aid, // Selection matrix for independent parameters
        const MatX& Ad, // Selection matrix for dependent and unidentifiable parameters
        const MatX& Kd, // Restructuring transformation matrix
        const double N,
        const double b_dim,
        const double d_dim,
        const double k, //hard redescender
        const double Weight_tor,
        const double Alpha_tor,
        const double O_sqrt_tol,
        const bool regroup,
        const bool includeOffset,
        Model &model_input
        // Eigen::VectorXi jtype_input
)
    :   pi_d_(pi_d), 
        Ginv_(Ginv), Aid_(Aid), Ad_(Ad), Kd_(Kd),
        N_(N), b_dim_(b_dim), d_dim_(d_dim), k_(k), Weight_tor_(Weight_tor), Alpha_tor_(Alpha_tor),O_sqrt_tol_(O_sqrt_tol),
        regroup_(regroup), includeOffset_(includeOffset)

{
    /* ----------------------------- set parameters ----------------------------- */
    // set the share point for the model 
    
    modelPtr_ = std::make_shared<Model>(model_input);
    nLinks_ = regroup_? modelPtr_->nv: 1;

    // Inertia parameters number
    p_ip_ =10 * nLinks_;
    p_ip_end_ = 10;

    // full parameters, b_full is regroup and p_full is not regroup
    const int fm_num = includeOffset_? 4:3;
    p_full_ = p_ip_ + fm_num * nLinks_;
    b_full_ = b_dim_ + fm_num * nLinks_;
    p_end_full_ = p_ip_end_ + fm_num * nLinks_;



    // file name 
    position="exciting-position.csv";
    velocity = "exciting-velocity.csv";
    acceleration ="exciting-acceleration.csv";
    solution ="exciting-solution.csv";

  


    if (regroup_){
        Wfull_ = MatX::Zero(modelPtr_->nv *N_, b_full_);
        alphanew_ = VecX::Ones(nLinks_);
    }else{
        Wfull_ = MatX::Zero(modelPtr_->nv *N_, p_end_full_);
        alphanew_ = VecX::Ones(nLinks_);
        std::cout <<"nLinks_*N_ " <<nLinks_*N_ << std::endl;
        std::cout <<"p_end_full_ " <<nLinks_*N_ << std::endl;

    }
    
    /* ------------  Initialize the full parameters from the model ------------ */

    // Structure [Inertia*nLinks_,Ia*nLinks_, firction *nLinks_]

    if (regroup_){
        X0_1_ = VecX::Zero(p_full_);
        for (int i = 0; i < nLinks_; i++) {
            const int pinocchio_joint_id = i + 1;
            MatX I = modelPtr_->inertias[pinocchio_joint_id].matrix();
            // intertia 
            X0_1_.segment<10>(10 * i)<< I(0,0), I(0,1), I(0,2),I(1,1),I(1,2),I(2,2),I(2,4),I(0,5),I(1,3),I(3,3);
            // Ia
            X0_1_.segment<1>(p_ip_+ i) << modelPtr_->rotorInertia(i);
            // friction 
            if (includeOffset_){
                X0_1_.segment(p_ip_ + nLinks_ + (fm_num-1) * i, fm_num-1) << modelPtr_->friction(i), modelPtr_->damping(i), 0.0;
            }else {
                X0_1_.segment(p_ip_ + nLinks_ + (fm_num-1) * i, fm_num-1) << modelPtr_->friction(i), modelPtr_->damping(i);
            }
        }
    }else{
        X0_1_ = VecX::Zero(p_end_full_);
        MatX I = modelPtr_->inertias[nLinks_].matrix();
        // std::cout << I << std::endl;
        
        // intertia 
        X0_1_.segment<10>(0)<< I(0,0), I(0,1), I(0,2),I(1,1),I(1,2),I(2,2),I(2,4),I(0,5),I(1,3),I(3,3);

        for (int i = 0; i < nLinks_; i++) {
            // const int pinocchio_joint_id = i + 1;
            // MatX I = modelPtr_->inertias[pinocchio_joint_id].matrix();
            // // intertia 
            // X0_1_.segment<10>(10 * i)<< I(0,0), I(0,1), I(0,2),I(1,1),I(1,2),I(2,2),I(2,4),I(0,5),I(1,3),I(3,3);
            // Ia
            X0_1_.segment<1>(p_ip_end_ + i) << modelPtr_->rotorInertia(i);
            // friction 
            if (includeOffset_){
                X0_1_.segment(p_ip_end_ + nLinks_ + (fm_num-1) * i, fm_num-1) << modelPtr_->friction(i), modelPtr_->damping(i), 0.0;
            }else {
                X0_1_.segment(p_ip_end_  + nLinks_ + (fm_num-1) * i, fm_num-1) << modelPtr_->friction(i), modelPtr_->damping(i);
            }
        }

    }


    /* ------------------------- lower and upper bounds ------------------------- */
    

    // Initialize lower and upper bounds, after regroup
    if (regroup_){
        VecX lowerB;
        VecX upperB;
        lb_ = VecX::Zero(b_full_);
        ub_ = VecX::Zero(b_full_);
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

        lb_.segment(b_dim_, fm_num * nLinks_).setConstant(0.0);
        ub_.segment(b_dim_, fm_num * nLinks_).setConstant(20.0);

    }else{
        lb_ = VecX::Zero(p_end_full_);
        ub_ = VecX::Zero(p_end_full_);

        lb_.segment(0,p_ip_end_) =X0_1_.segment(0,p_ip_end_) * 0.9;
        ub_.segment(0,p_ip_end_)  = X0_1_.segment(0,p_ip_end_) * 1.1;

        for (int j = 0; j < p_ip_end_; ++j) {
            if (lb_(j) >= ub_(j)) {
                std::swap(lb_(j), ub_(j));
            }
        }

        lb_.segment(p_ip_end_, fm_num * nLinks_).setConstant(0.0);
        ub_.segment(p_ip_end_, fm_num * nLinks_).setConstant(20.0);
    }
    
 



    
    /* ---------------------- load the optimized trajectory --------------------- */

    // using trajectories structure to save data
    // trajPtr_=std::make_shared<Trajectories>();
    // // later should move file name into sys_id input

    // importCSV(trajPtr_->q, "exciting-position.csv", 2000, nLinks_);
    // importCSV(trajPtr_->q_d, "exciting-velocity.csv", 2000, nLinks_);
    // importCSV(trajPtr_->q_dd, "exciting-acceleration.csv", 2000, nLinks_);


    // trajPtr_ = std::make_shared<FixedFrequencyFourierCurves>(10.0, 
    //                                                          50, 
    //                                                          model_input.nv, 
    //                                                          TimeDiscretization::Uniform, 
    //                                                          5,
    //                                                          2.0 * M_PI / 10);

     


    /* ------------------ process Observation matrix and Torque ----------------- */
    // ridPtr_ = std::make_shared<RegressorInverseDynamics>(model_input, trajPtr_,jtype_input);
    // do we need jtypep_input 

       // importCSV(trajPtr_->q, "exciting-position.csv", 2000, nLinks_);
    // importCSV(trajPtr_->q_d, "exciting-velocity.csv", 2000, nLinks_);
    // importCSV(trajPtr_->q_dd, "exciting-acceleration.csv", 2000, nLinks_);

    ridPtr_ = std::make_shared<RegressorInverseDynamics>(model_input,
                                                        position,
                                                        velocity,
                                                        acceleration);
    ridPtr_->computeWithDataImport(solution);


    // load z from solution 
    // z_.resize(13 * modelPtr_->nv);
    // Eigen::Array<VecX, 1, Eigen::Dynamic> data;
    // importCSV(data, "exciting-solution.csv",13 * modelPtr_->nv , 1);
    // // Array<VecX, 1, Eigen::Dynamic> to vector 
    // int index = 0;
    // for (int i = 0; i < data.size(); ++i) {
    //     for (int j = 0; j < data[i].size(); ++j) {
    //         z_(index++) = data[i](j);
    //     }
    // }
    // std::cout<<z<<std::endl;

    // ridPtr_->compute(z_, true);

}

//  may be should add in utils 
// void SysIDAlgFull::importCSV(Eigen::Array<VecX, 1, Eigen::Dynamic>& data,
void SysIDAlgFull::importCSV(Eigen::Array<VecX, 1, Eigen::Dynamic>& data,
               const std::string& filename,
               int expectedRows, int expectedCols) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file " + filename);
    }

    std::string line;
    std::vector<VecX> tempData;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        VecX rowData(expectedCols);
        double value;
        int col = 0;

        while (iss >> value && col < expectedCols) {
            rowData(col++) = value;
        }

        if (col != expectedCols) {
            throw std::runtime_error("Inconsistent number of columns in file");
        }

        tempData.push_back(rowData);
    }

    if (tempData.size() != expectedRows) {
        throw std::runtime_error("Unexpected number of rows in file");
    }

    data.resize(1, expectedRows);
    for (int i = 0; i < expectedRows; ++i) {
        data(0, i) = tempData[i];
    }
}

void SysIDAlgFull::runAlgorithm() {

    // get the regressor matrix(inertia) and torque
    // ridPtr_->compute(z_, true);
    // W_ip_ =ridPtr_->Y.leftCols(10*nLinks_);
    // T_ = W_ip_* X0_1_.head(10*nLinks_);

    // std::cout<<W_ip_<<std::endl;

    
    if (regroup_) {
        // ridPtr_->compute(z_, true);
        ridPtr_->computeWithDataImport(solution);
        W_ip_ =ridPtr_->Y.leftCols(10*nLinks_);
        T_ = W_ip_* X0_1_.head(10*nLinks_);
        regroupParameters();
    }else{
        ridPtr_->computeWithDataImport( solution);
        // ridPtr_->compute(z_, true);

        int startCol = (modelPtr_->nv - 1) * 10;
        W_ip_ =ridPtr_->Y.middleCols(startCol, 10);
        T_ = W_ip_* X0_1_.head(10);
    }

    //  initialize alpha
    VecX alphaold = regroup_ ? VecX::Zero(nLinks_) : VecX::Zero(1);
    VecX alphanew = alphanew_; 

    while ((alphaold -alphanew).norm() > Alpha_tor_) {
        // first identtifcation weight and base parameters
        runFirstOptimization();
        
        // second identifcation frictionexponent
        alphaold = alphanew_; 
        runSecondOptimization();
        // udpata alphanew
        alphanew = alphanew_;
        break;
    }

}


void SysIDAlgFull::regroupParameters() {

    // regroup w_ip
    Wb_ = W_ip_ * Aid_;

    // updata params
    VecX X0_ip = Aid_.transpose() * X0_1_.head(10 * nLinks_) + Kd_ * Ad_.transpose() * X0_1_.head(10 * nLinks_);
    VecX X0_fm = X0_1_.tail(X0_1_.size() - 10 * nLinks_);

    X0_1_.resize(b_full_);
    X0_1_.head(b_dim_) = X0_ip;
    X0_1_.tail(X0_fm.size()) = X0_fm;

    // VecX lower_B = Aid_.transpose() * lb_.head(10*nLinks_) 
}



void SysIDAlgFull::Wextension(const VecX& alpha, MatX& Wfm, MatX& Wm, MatX& Wf) {
    int wf_cols = includeOffset_? 3 * nLinks_ :2 * nLinks_;
    Wf = MatX::Zero(modelPtr_->nv * N_, wf_cols);
    Wm = MatX::Zero(modelPtr_->nv * N_, nLinks_);

    
    for (int i = 0; i < N_; ++i) {
        const VecX& q = ridPtr_->trajPtr_q(i);
        const VecX& q_d = ridPtr_->trajPtr_q_d(i);
        const VecX& q_dd = ridPtr_->trajPtr_q_dd(i);
        // std::cout<<"4.3.1"<<std::endl;
        for (int j =0; j< nLinks_ ; ++j){

            Wm(i * modelPtr_->nv + j, j) = q_dd(j); 
            // vel = q_d(j);
            
            if (includeOffset_) {
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
    fm_dim_ = Wfm.cols();
    Wfm << Wm, Wf;
}

void SysIDAlgFull::runFirstOptimization() {
    // Initialize Weighted matrix
    VecX Dnew;
    VecX Dold;

    Dnew = VecX::Ones(N_*nLinks_);
    Dold = VecX::Zero(N_*nLinks_);

    // if (regroup_){
    //     Dnew = VecX::Ones(N_*nLinks_);
    //     Dold = VecX::Zero(N_*nLinks_);
    // }else{
    //     Dnew = VecX::Ones(N_);
    //     Dold = VecX::Zero(N_);
    // }

    VecX Rh;

    MatX De;
    if(regroup_){
        De = Dnew.replicate(1,b_full_);
    }else{
        De = Dnew.replicate(1,p_end_full_);
    }


    MatX Wfm, Wm, Wf; 
    Wextension(alphanew_, Wfm, Wm, Wf);

    std::cout <<"Wfull_ " <<Wfull_.rows() << ","<< Wfull_.cols()<< std::endl;
    std::cout <<"W_ip_ " <<W_ip_.rows() << ","<< W_ip_.cols()<< std::endl;
    std::cout <<"Wfm " <<Wfm.rows() << ","<< Wfm.cols()<< std::endl;

    if(regroup_){
        Wfull_ << Wb_, Wfm;
    }else{
        Wfull_ << W_ip_ , Wfm;
    }


    // initialize covariance matrix
    MatX Oold ;
    MatX Onew ;

    Oold = MatX::Zero(nLinks_ , nLinks_ );
    Onew = MatX::Identity(nLinks_ , nLinks_ ); 

    // if (regroup_){
    //     Oold = MatX::Zero(nLinks_ , nLinks_ );
    //     Onew = MatX::Identity(nLinks_ , nLinks_ ); 
    // }else{
    //     Oold = MatX::Zero(1 , 1 );
    //     Onew = MatX::Identity(1 , 1 ); 
    // }


    while ((Dold - Dnew).norm() > Weight_tor_){
        // Cholesky 
        MatX O_sqrt = Onew.llt().matrixL();
        while ((Onew - Oold).norm() > O_sqrt_tol_){
            VecX Ts;
            MatX Ws;
            Ts = VecX::Zero(nLinks_*N_);
            Ws = MatX::Zero(nLinks_ * N_, Wfull_.cols());
            // MatX Wfs = MatX:Zero(nLinks_ * N_, Wfull.cols());
            for (int i =0; i < N_; ++i){
                Ts.segment(i*nLinks_, nLinks_) = O_sqrt.inverse()* T_.segment(i*nLinks_,nLinks_);
                Ws.block(i*nLinks_, 0, nLinks_, Wfull_.cols()) = O_sqrt.inverse() * Wfull_.block(i * nLinks_, 0, nLinks_, Wfull_.cols());
            }
            // if (regroup_){
            //     // normalized observation matrix and the observation response
            //     Ts = VecX::Zero(nLinks_*N_);
            //     Ws = MatX::Zero(nLinks_ * N_, Wfull_.cols());
            //     // MatX Wfs = MatX:Zero(nLinks_ * N_, Wfull.cols());
            //     for (int i =0; i < N_; ++i){
            //         Ts.segment(i*nLinks_, nLinks_) = O_sqrt.inverse()* T_.segment(i*nLinks_,nLinks_);
            //         Ws.block(i*nLinks_, 0, nLinks_, Wfull_.cols()) = O_sqrt.inverse() * Wfull_.block(i * nLinks_, 0, nLinks_, Wfull_.cols());
            //     }
            // }else{
            //     Ts = VecX::Zero(N_);
            //     Ws = MatX::Zero(N_, Wfull_.cols());
            //     // MatX Wfs = MatX:Zero(nLinks_ * N_, Wfull.cols());
            //     for (int i =0; i < N_; ++i){
            //         Ts.segment(i, 1) = O_sqrt.inverse()* T_.segment(i,1);
            //         Ws.block(i, 0, 1, Wfull_.cols()) = O_sqrt.inverse() * Wfull_.block(i , 0, 1, Wfull_.cols());
            //     }

            // }

           

            Th_ = Dnew.array() * Ts.array();
            Wh_ = De.array() * Ws.array(); //wfull

            //optimise funciton 
            SmartPtr<BaseParametersIdentification> mynlp = new BaseParametersIdentification();
            try {
                mynlp->BaseParametersIdentification::set_parameters(
                                    Ginv_,
                                    Wh_,
                                    Th_,
                                    pi_d_,
                                    X0_1_,
                                    b_dim_,
                                    fm_dim_,
                                    d_dim_,
                                    includeOffset_,
                                    *modelPtr_,
                                    lb_,
                                    ub_,
                                    p_ip_end_,
                                    regroup_);
                mynlp->constr_viol_tol = 1e-5;
                // mynlp->constr_viol_tol = 1;
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

            VecX X_sol = mynlp->solution;

            // updata initial 
            X0_1_ = X_sol;
            // calculate normlized residual
            // std::cout << "X0_1_: " << X0_1_.size() << std::endl;
            // std::cout << "Wh_: " << Wh_.rows()<< " ,"<< Wh_.cols() << std::endl;
            // std::cout << "Th_: " << Th_.size() << std::endl;
            Rh= Th_ - Wh_* X0_1_;
            // reshape Rh into n*m
            MatX Eh;
            if (regroup_){
                Eh = Eigen::Map<MatX>(Rh.data(), nLinks_, N_);
            }else{
                Eh = Eigen::Map<MatX>(Rh.data(), 1, N_);
            }
            // update covariance matrix

            Oold = Onew;
            if (regroup_){
                Onew = (Oold *Eh *Eh.transpose() * Oold)/(N_ - b_full_);
            }else{
                Onew = (Oold *Eh *Eh.transpose() * Oold)/(N_ - p_end_full_);
            }

            // update sqrt covariance matrix
            O_sqrt=  Onew.llt().matrixL();

            if (!O_sqrt.allFinite()) {
                break;
            }
            break;
        }
        // update weight 
        Dold = Dnew;
        // T clss hard redescender
        VecX Psi = (Rh.array().abs()<=k_).select(VecX::Ones(Rh.size()), VecX::Zero(Rh.size()));
        Dnew = Dold.array().min(Psi.array());

        if (regroup_) {
            De = Dnew.replicate(1, b_full_);
        } else {
            De = Dnew.replicate(1, p_end_full_);
        }
        break;
    }
}

void SysIDAlgFull::runSecondOptimization() {

    // MatX Wfm, Wm, Wf; // 只要一个？
    // Wentension(alphanew_, Wfm, Wm, Wf);
    // Wfull_ << Wb_, Wfm;

    // 
    // std::cout<<"wh"<<Wh_.cols()<<std::endl;
    int tailSize = includeOffset_ ? 4 * nLinks_ : 3 * nLinks_;

    VecX X0_2= VecX::Zero(10 * nLinks_ + tailSize);

    // if (regroup_){
    //     X0_2 = VecX::Zero(10 * nLinks_+ nLinks_ + tailSize);
    // }else{
    //     X0_2 = VecX::Zero(10  + nLinks_ + tailSize);
    // }

    int segmentSize = regroup_ ? b_dim_ + nLinks_ : 10 + nLinks_;
    X0_2.head(segmentSize) = X0_1_.head(segmentSize);

    // std::cout<<"wh"<< Wh_.rows()<<","<< Wh_.cols()<<std::endl;
    // std::cout<<"X0_2"<< X0_2.rows()<<std::endl;
    // std::cout<<"T_"<< T_.size()<<std::endl;
    // std::cout<<"T_"<< Th_.size()<<std::endl;
    VecX Fest = Th_ - Wh_ * X0_2;
    //  get last 4n params, only last 3 is useful
    VecX Xf = X0_1_.tail(4*nLinks_);
    Xf.head(nLinks_) = alphanew_;

    // if (regroup_){
    //     Xf.head(nLinks_) = alphanew_;
    // }else{
    //     Xf.head(1) = alphanew_;
    // } 
    // std::cout<<"Xf"<< Xf<<std::endl;

    VecX alphaold = VecX::Zero(nLinks_);
    VecX alphanew = alphanew_; 

    //optimise funciton 
    SmartPtr<FrictionParametersIdentification> mynlp_second = new FrictionParametersIdentification();
    try {
        mynlp_second->FrictionParametersIdentification::set_parameters(
                                            Xf,
                                            nLinks_,
                                            Fest,
                                            includeOffset_,
                                            N_,
                                            ridPtr_,
                                            regroup_);
                         
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

    std::cout <<"Xf "<<Xf<<std::endl;
    
    // update alphanew_
    if (regroup_){
        alphanew_ = Xf.head(nLinks_);
    }else{
        alphanew_ = Xf.head(1);
    }
    
    // update friction 
    if (includeOffset_) {
        if (regroup_) {
           X0_1_.segment(b_dim_ + nLinks_, 3 * nLinks_) = Xf.tail(3 * nLinks_);
        } else {
         X0_1_.segment(11 * nLinks_, 3 * nLinks_) = Xf.tail(3 * nLinks_);
        }
    } else {
        if (regroup_) {
            X0_1_.segment(b_dim_ + nLinks_, 2 * nLinks_) = Xf.tail(2 * nLinks_);
        } else {
           X0_1_.segment(10 + nLinks_, 2 * nLinks_) = Xf.tail(2 * nLinks_);
        }
    }


}

Eigen::VectorXd SysIDAlgFull::getfinalparam() const {
    return X0_1_;
}

Eigen::VectorXd SysIDAlgFull::getAlphaNew() const {
    return alphanew_;
}


} // namespace Kinova
} // namespace RAPTOR
