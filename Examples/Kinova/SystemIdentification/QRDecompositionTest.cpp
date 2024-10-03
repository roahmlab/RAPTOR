// main.cpp

#include "QRDecompositionSolver.h"
#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include "SysIDAlgFull.h"

int main() {

    int nq_nt = 100; 
    int p = 70;   

    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename, model);

    // model.gravity.linear()(2) = GRAVITY;
    model.friction.setZero();
    model.damping.setZero();
    model.rotorInertia.setZero(); 



    // Eigen::VectorXd X0_1_ = Eigen::VectorXd::Zero(10*model.nv);
    // for (int i = 0; i < model.nv; i++) {
    //     const int pinocchio_joint_id = i + 1;
    //     Eigen::MatrixXd I = modelPtr_->inertias[pinocchio_joint_id].matrix();
    //     // intertia 
    //     X0_1_.segment<10>(10 * i)<< I(0,0), I(0,1), I(0,2),I(1,1),I(1,2),I(2,2),I(2,4),I(0,5),I(1,3),I(3,3);
    //     // Ia
    //     // X0_1_.segment<1>(p_ip_+ i) << model.rotorInertia(i);
    //     // // friction 
    //     // if (includeOffset_){
    //     //     X0_1_.segment(p_ip_ + nLinks_ + (fm_num-1) * i, fm_num-1) << model.friction(i), model.damping(i), 0.0;
    //     // }else {
    //     //     X0_1_.segment(p_ip_ + nLinks_ + (fm_num-1) * i, fm_num-1) << model.friction(i), model.damping(i);
    //     // }
    // }


    RAPTOR::Kinova::QRDecompositionSolver qrSolver; 


    qrSolver.getData();


    // std::cout << "Aid matrix:\n" << qrSolver.Aid << "\n\n";
    // std::cout << "Ad matrix:\n" << qrSolver.Ad << "\n\n";
    // std::cout << "Kd matrix:\n" << qrSolver.Kd << "\n\n";
    // std::cout << "Beta vector:\n" << qrSolver.Beta << "\n\n";
    // std::cout << "Ginv matrix:\n" << qrSolver.Ginv << "\n\n";
    // std::cout << "Ginv RegroupMatrix:\n" << qrSolver.RegroupMatrix << "\n\n";
    // std::cout << "Dimension of independent parameters: " << qrSolver.dim_id << "\n";
    // std::cout << "Dimension of dependent parameters: " << qrSolver.dim_d << "\n";

    // Eigen::VectorXd phi1(qrSolver.Beta.size()+qrSolver.dim_d);
    // phi1 << qrSolver.Beta, qrSolver.pi_d;


    // Eigen::VectorXd phi2(qrSolver.Beta.size()+qrSolver.dim_d);
    // Eigen::VectorXd zero = Eigen::VectorXd::Zero(qrSolver.dim_d);
    // phi2 << qrSolver.Beta, zero;
    // std::cout << "groud : " <<qrSolver.InParam.transpose().head(14)<< "\n";
    // std::cout << "result1: " <<(qrSolver.Ginv *  phi1).transpose().head(14)<< "\n";
    // std::cout << "result2: " <<(qrSolver.Ginv *  phi2).transpose().head(14)<< "\n";
    double N =2000;
    double k = 3.0;
    double Weight_tor = 1e-3;
    double Alpha_tor= 1e-3;
    double O_sqrt_tol= 1e-3;


    // RAPTOR::Kinova::SysIDAlgFull sysIDSolver(
    //     qrSolver.pi_d, 
    //     qrSolver.Ginv,
    //     qrSolver.Aid,
    //     qrSolver.Ad,
    //     qrSolver.Kd,
    //     N,
    //     qrSolver.dim_id,
    //     qrSolver.dim_d,
    //     k,
    //     Weight_tor,
    //     Alpha_tor,
    //     O_sqrt_tol,
    //     false,
    //     true,
    //     model
    // );
    // sysIDSolver.runAlgorithm();
    // // sysIDSolver.runFirstOptimization();


 


    return 0;
}
