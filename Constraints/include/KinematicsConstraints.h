
#ifndef KINEMATICS_CONSTRAINTS_H
#define KINEMATICS_CONSTRAINTS_H

#include <unsupported/Eigen/MatrixFunctions>

#include "Utils.h"
#include "Constraints.h"
#include "Trajectories.h"
#include "ForwardKinematics.h"

#include <memory>

namespace IDTO {

inline Eigen::VectorXd flatRotationMatrix(const Eigen::Matrix3d& R) {
    Eigen::VectorXd x(9);
    x(0) = R(0, 0);
    x(1) = R(0, 1);
    x(2) = R(0, 2);
    x(3) = R(1, 0);
    x(4) = R(1, 1);
    x(5) = R(1, 2);
    x(6) = R(2, 0);
    x(7) = R(2, 1);
    x(8) = R(2, 2);
    return x;
}

class KinematicsConstraints : public Constraints {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    KinematicsConstraints() = default;

    KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                           const Model& model_input,
                           const Eigen::VectorXi& jtype_input,
                           const size_t joint_id_input,
                           const size_t time_id_input,
                           const Transform& desiredTransform_input,
                           const Transform endT_input = Transform());

    KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                           const Model& model_input,
                           const Eigen::VectorXi& jtype_input,
                           const size_t joint_id_input,
                           const size_t time_id_input,
                           const Vec3& desiredPosition_input,
                           const Transform endT_input = Transform());

    KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                           const Model& model_input,
                           const Eigen::VectorXi& jtype_input,
                           const size_t joint_id_input,
                           const size_t time_id_input,
                           const Mat3& desiredRotation_input,
                           const Transform endT_input = Transform());

    // Destructor
    ~KinematicsConstraints() = default;

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

    Vec3 desiredPosition;
    Mat3 desiredRotation;

    bool constrainPosition = false;
    bool constrainRotation = false;

    size_t joint_id = 0;

    size_t time_id = 0;

        // the transform matrix at the beginning and at the end
    Transform startT;
    Transform endT;

        // updated in compute()
    Transform jointT;
    MatX jointTJ;
    // Eigen::Array<MatX, 3, 1> jointTH;
    Eigen::Array<MatX, 12, 1> jointTH;

        // forward kinematics derivatives
    std::vector<Transform> dTdq;
    std::vector<std::vector<Transform>> ddTddq;
};

}; // namespace IDTO

#endif // KINEMATICS_CONSTRAINTS_H
