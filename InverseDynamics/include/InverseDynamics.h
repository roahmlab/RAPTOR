#ifndef INVERSEDYNAMICS_H
#define INVERSEDYNAMICS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "pinocchio/parsers/urdf.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"
#include "pinocchio/algorithm/crba.hpp"

#include <cmath>
#include <memory>
#include <cstdio>
#include <iostream>
#include <cstdlib>

namespace IDTO {

class InverseDynamics {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    InverseDynamics() = default;

    // Constructor
    InverseDynamics(const Model& model_input, int N_input);

    // Destructor
    ~InverseDynamics() = default;

    // class methods:
    virtual void compute(const Eigen::Array<VecX, 1, Eigen::Dynamic>& q, 
                         const Eigen::Array<VecX, 1, Eigen::Dynamic>& v, 
                         const Eigen::Array<VecX, 1, Eigen::Dynamic>& a,
                         bool compute_derivatives = true);

    // class members:
    std::unique_ptr<Model> modelPtr_;
    std::unique_ptr<Data> dataPtr_;

    int N = 0; // number of time instances in tspan

    MatX rnea_partial_dq;
    MatX rnea_partial_dv;
    MatX rnea_partial_da;

        // compute results are stored here
    Eigen::Array<VecX, 1, Eigen::Dynamic> tau;

        // compute results are stored here
    Eigen::Array<MatX, 1, Eigen::Dynamic> ptau_pq;
    Eigen::Array<MatX, 1, Eigen::Dynamic> ptau_pv;
    Eigen::Array<MatX, 1, Eigen::Dynamic> ptau_pa;
};

}; // namespace IDTO

#endif // INVERSEDYNAMICS_H
