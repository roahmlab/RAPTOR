#ifndef SYSIDALGCOMBINEALLCASES_H
#define SYSIDALGCOMBINEALLCASES_H

#include <vector>
#include <Eigen/Dense>
#include "Utils.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <cstdio>
#include <cstdlib>
#include "Optimizer.h"
#include "RegressorInverseDynamics.h"
#include "FixedFrequencyFourierCurves.h"
#include "QRDecompositionSolver.h"


// #include "JointLimits.h"
// #include "VelocityLimits.h"
// #include "TorqueLimits.h"
// #include "KinovaCustomizedConstraints.h"

namespace RAPTOR {
namespace Kinova {

class SysIDAlgCombineAllCases{
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Model = pinocchio::Model;

    SysIDAlgCombineAllCases(
        const int N,
        const int Alg_case,
        const double error_threshold, //hard redescender
        const double Weight_tol,
        const double Alpha_tol,
        const double O_sqrt_tol,
        const bool include_friction_offset,
        std::vector<std::string> files,
        Model &model_input
    );

    void runAlgorithm();
    VecX getfinalparam() const;
    VecX getAlphaNew() const;



private:
    // input
    const int N_;
    const int Alg_case_;
    const double error_threshold_;
    const double Weight_tol_;
    const double Alpha_tol_;
    const double O_sqrt_tol_;
    const bool include_friction_offset_;

    // file name 
    std::string position;
    std::string velocity;
    std::string acceleration;
    std::string solution;

    // regroup
    VecX pi_d_;
    MatX Ginv_;
    MatX Aid_;
    MatX Ad_;
    MatX Kd_;
    int b_dim_;
    int d_dim_;
    std::shared_ptr<QRDecompositionSolver> regroupPtr_;

    // set parameters
    int fm_num_;
    int nLinks_; 
    int p_ip_;
    int p_full_;
    int b_full_;
    MatX Wfull_; 
    VecX alphanew_;

    //  model, boundary,Observation matrix and Torque
    VecX X0_1_;
    VecX lb_;
    VecX ub_;
    VecX T_;
    MatX W_ip_;
    MatX Wb_;
    std::shared_ptr<Model> modelPtr_;
    std::shared_ptr<RegressorInverseDynamics> ridPtr_;

    VecX Th_;
    MatX Wh_;

   
    void regroupParameters();
    void runFirstOptimization();
    void runSecondOptimization(VecX& alphaold, VecX& alphanew);
    void Wextension(const VecX& alpha, MatX& Wfm, MatX& Wm, MatX& Wf);

};

} // namespace Kinova
} // namespace RAPTOR



#endif //SYSIDALGCOMBINEALLCASES_H