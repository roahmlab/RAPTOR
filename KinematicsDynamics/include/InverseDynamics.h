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
#include "pinocchio/algorithm/rnea-second-order-derivatives.hpp"
#include "pinocchio/algorithm/crba.hpp"

#include "Trajectories.h"

#include <cmath>
#include <memory>
#include <cstdio>
#include <iostream>
#include <cstdlib>

namespace RAPTOR {

class InverseDynamics {
public:
    using Model = pinocchio::ModelTpl<double>;
    using Data = pinocchio::DataTpl<double>;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;
    using Ten3 = pinocchio::DataTpl<double>::Tensor3x;

    // Constructor
    InverseDynamics() = default;

    // Constructor
    InverseDynamics(const Model& model_input, 
                    const std::shared_ptr<Trajectories>& trajPtr_input);

    // Destructor
    ~InverseDynamics() = default;

    // class methods:
    virtual void compute(const VecX& z,
                         bool compute_derivatives = true,
                         bool compute_hessian = false);

    MatX chipFromTensor3x(const Ten3& tensor3x, 
                          const Eigen::Index offset, 
                          const Eigen::Index dim);

        // determine if the constraints are computed before and save the current decision variable
    bool is_computed(const VecX& z, 
                     bool compute_derivatives, 
                     bool compute_hessian);

    // class members:
    std::shared_ptr<Model> modelPtr_ = nullptr;
    std::shared_ptr<Data> dataPtr_ = nullptr;

    std::shared_ptr<Trajectories> trajPtr_ = nullptr;

    int N = 0; // number of time instances in tspan

        // the decision variable that was evaluated last time
    VecX current_z;
    bool if_compute_derivatives = false;
    bool if_compute_hessian = false;

    MatX prnea_pq;
    MatX prnea_pv;
    MatX prnea_pa;

    Ten3 ptau2_pq;
    Ten3 ptau2_pv;
    Ten3 ptau2_pqpv;
    Ten3 ptau2_papq;

        // compute results are stored here
    Eigen::Array<VecX, 1, Eigen::Dynamic> tau;
    Eigen::Array<MatX, 1, Eigen::Dynamic> ptau_pz;
    Eigen::Array<MatX, Eigen::Dynamic, Eigen::Dynamic> ptau_pz_pz;
};

}; // namespace RAPTOR

#endif // INVERSEDYNAMICS_H
