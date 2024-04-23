#ifndef CUSTOMIZED_INVERSEDYNAMICS_H
#define CUSTOMIZED_INVERSEDYNAMICS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "pinocchio/parsers/urdf.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"

#include "InverseDynamics.h"
#include "Spatial.h"
#include "Trajectories.h"

#include <cmath>
#include <memory>
#include <cstdio>
#include <iostream>
#include <cstdlib>

namespace IDTO {

// rewrite the inverse dynamics using the original Roy Featherstone's algorithm
// so that we can get the force and the gradient of the force on the fixed joint
class CustomizedInverseDynamics : public InverseDynamics {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Vec6 = Vector6d;
    using Mat6 = Matrix6d;

    // Constructor
    CustomizedInverseDynamics() = default;

    // Constructor
    CustomizedInverseDynamics(const Model& model_input, 
                              const Eigen::VectorXi& jtype_input,
                              const std::shared_ptr<Trajectories>& trajPtr_input);

    // Destructor
    ~CustomizedInverseDynamics() = default;

    // class methods:
    virtual void compute(const VecX& z,
                         bool compute_derivatives = true); 

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
    Eigen::Array<Vec6, 1, Eigen::Dynamic> f;
    Eigen::Array<MatX, 1, Eigen::Dynamic> pf_pz;

        // compute results are stored here
    Eigen::Array<Vec6, 1, Eigen::Dynamic> lambda;
    Eigen::Array<MatX, 1, Eigen::Dynamic> plambda_pz;
};

}; // namespace IDTO

#endif // CUSTOMIZED_INVERSEDYNAMICS_H
