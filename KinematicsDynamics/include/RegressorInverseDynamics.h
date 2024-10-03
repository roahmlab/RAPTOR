#ifndef REGRESSOR_INVERSEDYNAMICS_H
#define REGRESSOR_INVERSEDYNAMICS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include <pinocchio/algorithm/regressor.hpp>


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
    using Vec6 = Vector6d;
    using Mat6 = Matrix6d;

    // Constructor
    RegressorInverseDynamics() = default;

    // Constructor
    RegressorInverseDynamics(const Model& model_input, 
                             const std::shared_ptr<Trajectories>& trajPtr_input,
                             Eigen::VectorXi jtype_input = Eigen::VectorXi(0));

    RegressorInverseDynamics(const Model& model_input,
                            const std::string& position,
                            const std::string& velocity,
                            const std::string& acceleration, 
                            Eigen::VectorXi jtype_input= Eigen::VectorXi(0));

    // Destructor
    ~RegressorInverseDynamics() = default;

    // class methods:
    virtual void compute(const VecX& z,
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override; 

    virtual void computeWithDataImport(const std::string& solution);
                                           
                         
    // class members:
    Eigen::VectorXi jtype;
    Eigen::Array<Mat6, 1, Eigen::Dynamic> Xtree;
    Eigen::Array<Mat6, 1, Eigen::Dynamic> I;
    Vec6 a_grav;
    
    Eigen::Array<Mat6, 1, Eigen::Dynamic> Xup;
    Eigen::Array<Mat6, 1, Eigen::Dynamic> dXupdq;
    Eigen::Array<Vec6, 1, Eigen::Dynamic> S;
    Eigen::Array<Vec6, 1, Eigen::Dynamic> v;
    Eigen::Array<MatX, 1, Eigen::Dynamic> pv_pz;
    Eigen::Array<Vec6, 1, Eigen::Dynamic> a;
    Eigen::Array<MatX, 1, Eigen::Dynamic> pa_pz;
    
    VecX phi;

    int numInertialParams = 0;
    int numParams = 0;

    MatX Y;
    Eigen::Array<MatX, 1, Eigen::Dynamic> pY_pz;

    MatX Yfull;
    Eigen::Array<MatX, 1, Eigen::Dynamic> pYfull_pz;

    MatX Ycurrent;

    // load data
    Eigen::Array<VecX, 1, Eigen::Dynamic> trajPtr_q;
    Eigen::Array<VecX, 1, Eigen::Dynamic> trajPtr_q_d;
    Eigen::Array<VecX, 1, Eigen::Dynamic> trajPtr_q_dd;
    // MatX trajPtr_q;
    // MatX trajPtr_q_d;
    // MatX trajPtr_q_dd;
};

}; // namespace RAPTOR

#endif // REGRESSOR_INVERSEDYNAMICS_H