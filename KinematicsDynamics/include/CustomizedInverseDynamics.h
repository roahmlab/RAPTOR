#ifndef CUSTOMIZED_INVERSEDYNAMICS_H
#define CUSTOMIZED_INVERSEDYNAMICS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include "InverseDynamics.h"
#include "Spatial.h"
#include "Transform.h"
#include "Trajectories.h"

#include <cmath>
#include <memory>
#include <cstdio>
#include <iostream>
#include <cstdlib>

namespace RAPTOR {

// rewrite the inverse dynamics using the original Roy Featherstone's algorithm
// so that we can get the force and the gradient of the force on the fixed joint
class CustomizedInverseDynamics : public InverseDynamics {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Vec6 = Eigen::Vector<double, 6>;
    using Mat6 = Eigen::Matrix<double, 6, 6>;

    // Constructor
    CustomizedInverseDynamics() = default;

    // Constructor
    CustomizedInverseDynamics(const Model& model_input, 
                              const std::shared_ptr<Trajectories>& trajPtr_input,
                              Eigen::VectorXi jtype_input = Eigen::VectorXi(0));

    // Destructor
    ~CustomizedInverseDynamics() = default;

    // class methods:
    VecX get_full_joints(const VecX& q) const;

    MatX get_full_joints_derivative(const MatX& q) const;

    virtual void compute(const VecX& z,
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override; 

    // class members:
    Eigen::VectorXi jtype;
    std::vector<int> active_joints;

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
    Eigen::Array<Vec6, 1, Eigen::Dynamic> f;
    Eigen::Array<MatX, 1, Eigen::Dynamic> pf_pz;

        // compute results are stored here
    Eigen::Array<Vec6, 1, Eigen::Dynamic> lambda;
    Eigen::Array<MatX, 1, Eigen::Dynamic> plambda_pz;
};

}; // namespace RAPTOR

#endif // CUSTOMIZED_INVERSEDYNAMICS_H
