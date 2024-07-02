#ifndef KINOVA_CUSTOMIZED_CONSTRAINTS_H
#define KINOVA_CUSTOMIZED_CONSTRAINTS_H

#include "Constraints.h"
#include "Trajectories.h"
#include "ForwardKinematics.h"
#include "BoxCollisionAvoidance.h"

#include <memory>

namespace IDTO {
namespace Kinova {

/*
This class implements an example of customized constraints for Kinova.
It uses spheres to represent forward occupancy of the robot, 
computes the distance between the robot and the customized obstacles,
and make sure the distance are larger than 0 to achieve collision avoidance.
*/

#define NUM_SPHERES 14

// based on which joint we compute the positions of the spheres
const int sphere_joint_id[NUM_SPHERES] = {2, 2, 2, 2, 2, 2, 2, // spheres on link 2
                                          4, 4, 4, 4, 4,       // spheres on link 4
                                          6, 6};               // spheres on link 6

// the translation axis of the spheres to cover the links (they are all extended along the y axis)
const int sphere_offset_axis[NUM_SPHERES] = {5, 5, 5, 5, 5, 5, 5, // spheres on link 2
                                             5, 5, 5, 5, 5,       // spheres on link 4
                                             5, 5};               // spheres on link 6

// the translation offset of the spheres to cover the links                                            
const double sphere_offset[NUM_SPHERES] = {-0.05, -0.11, -0.17, -0.24, -0.30, -0.36, -0.43, // spheres on link 2
                                           -0.07, -0.14, -0.21, -0.28, -0.32,               // spheres on link 4
                                           -0.07, -0.14};                                   // spheres on link 6                                        

const double sphere_radius[NUM_SPHERES] = {0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.08, // spheres on link 2
                                           0.07, 0.06, 0.06, 0.06, 0.07,             // spheres on link 4
                                           0.05, 0.05};                              // spheres on link 6

class KinovaCustomizedConstraints : public Constraints {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    KinovaCustomizedConstraints() = default;

    // Constructor
    KinovaCustomizedConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                const Model& model_input,
                                const Eigen::VectorXi& jtype_input,
                                const std::vector<Vec3>& boxCenters_input,
                                const std::vector<Vec3>& boxOrientation_input,
                                const std::vector<Vec3>& boxSize_input);

    // Destructor
    ~KinovaCustomizedConstraints() = default;

    // class methods:
        // compute constraints
    void compute(const VecX& z, 
                 bool compute_derivatives = true,
                 bool compute_hessian = false) override;

        // compute constraints lower bounds and upper bounds
    void compute_bounds() override;

        // print violation information
    void print_violation_info() override;

    // class variables:
    std::shared_ptr<Trajectories>& trajPtr_;

    std::shared_ptr<BoxCollisionAvoidance> collisionAvoidancePtr_;

    std::unique_ptr<Model> modelPtr_;

    std::unique_ptr<ForwardKinematicsSolver> fkhofPtr_;

        // jtype copy
    Eigen::VectorXi jtype;

        // the transform matrix at the beginning and at the end
    Transform startT;
    Transform endT;

        // updated in compute()
    Transform jointT;
    MatX jointTJ;
    MatX pq_pz;

        // forward kinematics derivatives
    std::vector<Transform> dTdq;
};

}; // namespace Kinova
}; // namespace IDTO

#endif // KINOVA_CUSTOMIZED_CONSTRAINTS_H
