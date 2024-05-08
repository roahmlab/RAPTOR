
#ifndef ENDEFFECTORCONSTRAINTS_H
#define ENDEFFECTORCONSTRAINTS_H

#include "Constraints.h"
#include "Trajectories.h"
#include "ForwardKinematics.h"

#include <memory>

namespace IDTO {

class EndEffectorConstraints : public Constraints {
public:
    using Model = pinocchio::Model;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    EndEffectorConstraints() = default;

    // Constructor
    EndEffectorConstraints(const Model& model_input,
                           const Eigen::VectorXi& jtype_input,
                           const Transform& endT_input,
                           const std::string joint_name_input,
                           std::shared_ptr<Trajectories>& trajPtr_input,
                           const VecX& desiredEndEffectorPos_input);

    // Destructor
    ~EndEffectorConstraints() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() override;

    // class members:
    std::shared_ptr<Trajectories>& trajPtr_;

    std::unique_ptr<Model> modelPtr_;
    Eigen::VectorXi jtype;

    std::unique_ptr<ForwardKinematicsHighOrderDerivative> fkhofPtr_;

    VecX desiredEndEffectorPos;

    Model::JointIndex joint_id = 0;

        // the transform matrix at the beginning and at the end
    Transform startT;
    Transform endT;

        // updated in compute()
    Transform jointT;
    MatX jointTJ;
    Eigen::Array<MatX, 6, 1> jointTH;

        // forward kinematics derivatives
    std::vector<Transform> dTdq;
    std::vector<std::vector<Transform>> ddTddq;
};

}; // namespace IDTO

#endif // ENDEFFECTORCONSTRAINTS_H
