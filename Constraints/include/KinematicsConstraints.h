#ifndef KINEMATICS_CONSTRAINTS_H
#define KINEMATICS_CONSTRAINTS_H

#include "Constraints.h"
#include "Trajectories.h"
#include "DynamicsConstraints.h"
#include "ForwardKinematics.h"

namespace IDTO {

class KinematicsConstraints : public Constraints {
public:
    using Model = pinocchio::Model;
    using Data = pinocchio::Data;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    KinematicsConstraints() = default;

    // Constructor
    KinematicsConstraints(const Model& model_input,
                          const Eigen::VectorXi& jtype_input,
                          std::shared_ptr<Trajectories>& trajPtr_input,
                          const std::string joint_name_input,
                          const MatX& lowerLimits_input,
                          const MatX& upperLimits_input,
                          const Transform startT_input = Transform(),
                          const Transform endT_input = Transform(),
                          std::shared_ptr<DynamicsConstraints> dcPtr_input = nullptr);

    // Destructor
    ~KinematicsConstraints() = default;

    // class methods:
        // compute constraints
    void compute(const VecX& z, bool compute_derivatives = true) override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() override;

    // class variables:
    std::shared_ptr<Trajectories> trajPtr_;
    std::shared_ptr<DynamicsConstraints> dcPtr_;

    // class members:
    std::unique_ptr<Model> modelPtr_;

    std::unique_ptr<ForwardKinematicsHighOrderDerivative> fkhofPtr_;

        // jtype copy
    Eigen::VectorXi jtype;

        // the joint index of the joint we want to constrain
    Model::JointIndex joint_id = 0;

        // the transform matrix at the beginning and at the end
    Transform startT;
    Transform endT;

        // updated in compute()
    Transform jointT;
    MatX jointTJ;
    MatX pq_pz;

        // forward kinematics derivatives
    std::vector<Transform> dTdq;

    MatX lowerLimits;
    MatX upperLimits;
};

}; // namespace IDTO

#endif // KINEMATICS_CONSTRAINTS_H
