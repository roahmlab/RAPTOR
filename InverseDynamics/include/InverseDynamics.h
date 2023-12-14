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

#include "Trajectories.h"

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
    InverseDynamics(const Model& model_input, 
                    std::shared_ptr<Trajectories>& trajPtr_input);

    // Destructor
    ~InverseDynamics() = default;

    // class methods:
    virtual void compute(const VecX& z,
                         bool compute_derivatives = true);

        // determine if the constraints are computed before and save the current decision variable
    bool is_computed(const VecX& z, bool compute_derivatives) {
        if (z.size() != current_z.size()) {
            current_z = z;
            if_compute_derivatives = compute_derivatives;
            return false;
        }

        if (compute_derivatives != if_compute_derivatives) {
            current_z = z;
            if_compute_derivatives = compute_derivatives;
            return false;
        }

        bool isSame = true;

        for (int i = 0; i < z.size(); i++) {
            if (z(i) != current_z(i)) {
                isSame = false;
                break;
            }
        }

        current_z = z;  
        if_compute_derivatives = compute_derivatives;
        return isSame;
    }                     

    // class members:
    std::unique_ptr<Model> modelPtr_;
    std::unique_ptr<Data> dataPtr_;

    std::shared_ptr<Trajectories> trajPtr_;

    int N = 0; // number of time instances in tspan

        // the decision variable that was evaluated last time
    VecX current_z;
    bool if_compute_derivatives = false;

    MatX prnea_pq;
    MatX prnea_pv;
    MatX prnea_pa;

        // compute results are stored here
    Eigen::Array<VecX, 1, Eigen::Dynamic> tau;

        // compute results are stored here
    Eigen::Array<MatX, 1, Eigen::Dynamic> ptau_pz;
};

}; // namespace IDTO

#endif // INVERSEDYNAMICS_H
