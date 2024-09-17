
#ifndef KINEMATICS_CONSTRAINTS_H
#define KINEMATICS_CONSTRAINTS_H

#include "Constraints.h"
#include "Trajectories.h"
#include "ForwardKinematics.h"
#include "LieSpaceResidual.h"

namespace RAPTOR {

class KinematicsConstraints : public Constraints {
public:
    using Model = pinocchio::ModelTpl<float>;
    using Vec3 = Eigen::Vector3f;
    using Mat3 = Eigen::Matrix3f;
    using VecX = Eigen::VectorXf;
    using MatX = Eigen::MatrixXf;

    // Constructor
    KinematicsConstraints() = default;

    KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                          const Model* model_input,
                          const size_t joint_id_input,
                          const size_t time_id_input,
                          const Transform& desiredTransform_input,
                          const Transform endT_input = Transform(),
                          Eigen::VectorXi jtype_input = Eigen::VectorXi(0));

    KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                          const Model* model_input,
                          const size_t joint_id_input,
                          const size_t time_id_input,
                          const Vec3& desiredPosition_input,
                          const Transform endT_input = Transform(),
                          Eigen::VectorXi jtype_input = Eigen::VectorXi(0));

    KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                          const Model* model_input,
                          const size_t joint_id_input,
                          const size_t time_id_input,
                          const Mat3& desiredRotation_input,
                          const Transform endT_input = Transform(),
                          Eigen::VectorXi jtype_input = Eigen::VectorXi(0));

    // Destructor
    ~KinematicsConstraints() = default;

    // class methods:
        // compute constraints
    virtual void compute(const VecX& z, 
                         bool compute_derivatives = true,
                         bool compute_hessian = false) override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() override;

    void print_violation_info() override;

    // class members:
    std::shared_ptr<Trajectories>& trajPtr_;

    const Model* modelPtr_ = nullptr;

    std::unique_ptr<ForwardKinematicsSolver> fkPtr_;

    Vec3 desiredPosition;
    Mat3 desiredRotation;

    bool constrainPosition = false;
    bool constrainRotation = false;

    size_t joint_id = 0;
    size_t time_id = 0;

        // the transform matrix at the beginning and at the end
    Transform startT;
    Transform endT;
};

}; // namespace RAPTOR

#endif // KINEMATICS_CONSTRAINTS_H
