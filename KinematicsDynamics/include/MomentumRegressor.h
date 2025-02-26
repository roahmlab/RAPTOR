#ifndef MOMENTUM_REGRESSOR_H
#define MOMENTUM_REGRESSOR_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/regressor.hpp"

#include "RegressorInverseDynamics.h"

#include <cmath>
#include <iostream> 
#include <memory>
#include <cstdio>
#include <cstdlib>

namespace RAPTOR {

// Compute system momentum using H(q) * qdot = Y * phi,
// where Y is the n x 10*n regressor matrix and 
// phi is the vector of 10*n dynamic parameters (inertia, com, mass for each of the link).
// phi is constant and directly loaded from the robot. 
// The gradient of Y will also be computed.
// Note that although MomentumRegressor is subclass of InverseDynamics, 
// it does not compute torque but the system momentum, which is stored in tau.
class MomentumRegressor : public RegressorInverseDynamics {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using Vec6 = Eigen::Vector<double, 6>;
    using Mat6 = Eigen::Matrix<double, 6, 6>;

    // Constructor
    MomentumRegressor() = default;

    // Constructor
    MomentumRegressor(const Model& model_input, 
                      const std::shared_ptr<Trajectories>& trajPtr_input,
                      Eigen::VectorXi jtype_input = Eigen::VectorXi(0)) :
        RegressorInverseDynamics(model_input, trajPtr_input, false, jtype_input) {
        Y_CTv.resize(N * modelPtr_->nv, numParams);
        Y_CTv.setZero();
        pY_CTv_pz.resize(trajPtr_->varLength);
        for (int i = 0; i < trajPtr_->varLength; i++) {
            pY_CTv_pz(i).resize(N * modelPtr_->nv, numParams);
            pY_CTv_pz(i).setZero();
        }
    };

    // Destructor
    ~MomentumRegressor() = default;

    // class methods:
    virtual void compute(const VecX& z,
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

    // class members:

    // this has been defined in RegressorInverseDynamics 
    // and stores the regressor for system momentum H(q) * v
    // MatX Y;
    // Eigen::Array<MatX, 1, Eigen::Dynamic> pY_pz;

    // this stores the regressor for C^T(q, v) * v
    // which is needed to compute the time derivative of the system momentum
    MatX Y_CTv;
    Eigen::Array<MatX, 1, Eigen::Dynamic> pY_CTv_pz;
};

}; // namespace RAPTOR

#endif // MOMENTUM_REGRESSOR_H