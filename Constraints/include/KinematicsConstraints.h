#ifndef KINEMATICS_CONSTRAINTS_H
#define KINEMATICS_CONSTRAINTS_H

#include "Constraints.h"
#include "Trajectories.h"
#include "ForwardKinematics.h"

namespace IDTO {

class KinematicsConstraints : public Constraints {
public:
    // Constructor
    KinematicsConstraints() = default;

    // Constructor
    KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                       const Eigen::VectorXi& jtype_input);

    // Destructor
    ~KinematicsConstraints() = default;

    // class methods:
        // compute constraints
    void compute(const VecX& z, bool compute_derivatives = true) override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() override;

    // class variables:
    std::shared_ptr<Trajectories> trajPtr_;

    // class members:
    std::unique_ptr<ForwardKinematicsHighOrderDerivative> fkhofPtr_;

        // jtype copy
    Eigen::VectorXi jtype;

        // for contact constraints (forward kinematics mainly)
    Model::JointIndex contact_joint_id = 0;

        // forward kinematics derivatives
    std::vector<Transform> dTdq;
    
};

}; // namespace IDTO

#endif // KINEMATICS_CONSTRAINTS_H
