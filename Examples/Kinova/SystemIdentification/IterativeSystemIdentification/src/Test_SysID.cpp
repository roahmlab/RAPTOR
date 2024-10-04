// main.cpp

#include "QRDecompositionSolver.h"
#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include "SysIDAlgCombineAllCases.h"

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


    int N =2000;
    double error_threshold = 3.0;
    double Weight_tol = 1;
    double Alpha_tol= 1e-3;
    double O_sqrt_tol= 1e-3;
    int Alg_case = 1; // 0 is full 1 is end 2 is f
    bool include_friction_offset = true;

    // file name 
    std::string position= "exciting-position.csv";
    std::string velocity = "exciting-velocity.csv";
    std::string acceleration ="exciting-acceleration.csv";
    std::string solution ="exciting-solution.csv";
    std::vector<std::string> files = {position, velocity, acceleration, solution};

    RAPTOR::Kinova::SysIDAlgCombineAllCases sysIDSolver(
        N,
        Alg_case,
        error_threshold,
        Weight_tol,
        Alpha_tol,
        O_sqrt_tol,
        include_friction_offset,
        files,
        model
    );
    sysIDSolver.runAlgorithm();
    // sysIDSolver.runFirstOptimization();

    return 0;
}
