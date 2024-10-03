#ifndef SYSIDALGFULL_H
#define SYSIDALGFULL_H

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


// #include "JointLimits.h"
// #include "VelocityLimits.h"
// #include "TorqueLimits.h"
// #include "KinovaCustomizedConstraints.h"

namespace RAPTOR {
namespace Kinova {

class SysIDAlgFull{
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Model = pinocchio::Model;

    SysIDAlgFull(
        // const VecX& X0_1, // Initial parameter estimation, p_full
        // const VecX& T, // Torque vector. n * m
        const VecX& pi_d,
        // const MatX& W_ip, // Observation matrix of the inertia parameters
        const MatX& Ginv, // The inverse of the bijective mapping that converts base inertia parameters into full inertia parameters
        const MatX& Aid, // Selection matrix for independent parameters
        const MatX& Ad, // Selection matrix for dependent and unidentifiable parameters
        const MatX& Kd, // Restructuring transformation matrix
        const double N,
        const double b_dim,
        const double d_dim,
        const double k, // hard redescender
        const double Weight_tor,
        const double Alpha_tor,
        const double O_sqrt_tol,
        const bool regroup,
        const bool includeOffset,
        Model &model_input
    );

    void runAlgorithm();
    VecX getfinalparam() const;
    VecX getAlphaNew() const;
    void runFirstOptimization();


private:

    int nLinks_;
    int N_;
    int b_dim_;
    int d_dim_;
    int fm_dim_;
    VecX lb_;
    VecX ub_;
    // VecX lba_;
    // VecX uba_;
    VecX X0_1_;
    VecX T_;
    MatX W_ip_;
    MatX Wb_;
    MatX Ginv_;
    MatX Aid_;
    MatX Ad_;
    MatX Kd_;
    const MatX pi_d_;

    // 
    double k_;
    double tol_;
    double Weight_tor_; 
    double O_sqrt_tol_; 
    double Alpha_tor_;
    bool regroup_;
    bool includeOffset_;


    // mid
    int p_ip_;
    int p_full_;
    int b_full_;
    int p_end_full_;
    int p_ip_end_;
    MatX Wfull_; 
    VecX alphanew_;
    VecX Th_;
    MatX Wh_;
    VecX z_;

    // file name 
    std::string position;
    std::string velocity;
    std::string acceleration;
    std::string solution;

   
    void regroupParameters();
    // void runFirstOptimization();
    void runSecondOptimization();
    void importCSV(Eigen::Array<VecX, 1, Eigen::Dynamic>& data,
                  const std::string& filename,
                  int N_, int nLinks_);
    void Wextension(const VecX& alpha, MatX& Wfm, MatX& Wm, MatX& Wf);

    std::shared_ptr<RegressorInverseDynamics> ridPtr_;
    std::shared_ptr<Trajectories> trajPtr_;
    std::shared_ptr<Model> modelPtr_;
};

} // namespace Kinova
} // namespace RAPTOR



#endif //SYSIDALGFULL_H