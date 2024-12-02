#ifndef REGRESSOR_INVERSEDYNAMICS_H
#define REGRESSOR_INVERSEDYNAMICS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/regressor.hpp"

#include "CustomizedInverseDynamics.h"
#include "Spatial.h"
#include "Trajectories.h"

#include <cmath>
#include <iostream> 
#include <memory>
#include <cstdio>
#include <cstdlib>

namespace RAPTOR {

// Compute inverse dynamics using tau = Y * phi,
// where Y is the n x 10*n regressor matrix and 
// phi is the vector of 10*n dynamic parameters (inertia, com, mass for each of the link).
// phi is constant and directly loaded from the robot. 
// The gradient of Y will also be computed.
class RegressorInverseDynamics : public InverseDynamics {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using Vec6 = Eigen::Vector<double, 6>;
    using Mat6 = Eigen::Matrix<double, 6, 6>;
    using MatRegressor = Eigen::Matrix<double, 6, 10>;

    // Constructor
    RegressorInverseDynamics() = default;

    // Constructor
    RegressorInverseDynamics(const Model& model_input, 
                             const std::shared_ptr<Trajectories>& trajPtr_input,
                             const bool include_motor_dynamics = true,
                             Eigen::VectorXi jtype_input = Eigen::VectorXi(0));

    // Destructor
    ~RegressorInverseDynamics() = default;

    // class methods:
    virtual void compute(const VecX& z,
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;
                                           
    // class members:
    Eigen::VectorXi jtype;
    Eigen::Array<Mat6, 1, Eigen::Dynamic> Xtree;
    Vec6 a_grav;
    
    VecX phi;

    int numInertialParams = 0;
    int numParams = 0;

    MatX Y;
    Eigen::Array<MatX, 1, Eigen::Dynamic> pY_pz;

    int NB = 0;
};

}; // namespace RAPTOR

#endif // REGRESSOR_INVERSEDYNAMICS_H