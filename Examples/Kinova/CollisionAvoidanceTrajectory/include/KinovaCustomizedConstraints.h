#ifndef KINOVA_CUSTOMIZED_CONSTRAINTS_H
#define KINOVA_CUSTOMIZED_CONSTRAINTS_H

#include "Constraints.h"
#include "Trajectories.h"
#include "ForwardKinematics.h"
#include "BoxCollisionAvoidance.h"

#include <memory>

namespace RAPTOR {
namespace Kinova {

/*
This class implements an example of customized constraints for Kinova.
It uses spheres to represent forward occupancy of the robot, 
computes the distance between the robot and the customized obstacles,
and make sure the distance are larger than 0 to achieve collision avoidance.
*/

constexpr int NUM_SPHERES = 17;

// based on which joint we compute the positions of the spheres
const int SPHERE_JOINT_ID[NUM_SPHERES] = {2, 2, 2, 2, 2, 2, 2, // spheres on link 2
                                          4, 4, 4, 4, 4,       // spheres on link 4
                                          6, 6,                // spheres on link 6
                                          7, 7, 7};            // spheres on camera

// the translation offset of the spheres to cover the links                                            
const double SPHERE_OFFSET[NUM_SPHERES][3] = {
    {0.0, -0.05, -0.01},
    {0.0, -0.11, -0.01},
    {0.0, -0.17, -0.01},
    {0.0, -0.24, -0.01},
    {0.0, -0.30, -0.01},
    {0.0, -0.36, -0.01},
    {0.0, -0.43, -0.01}, // sphere on link 2
    {0.0, -0.07, 0.0},
    {0.0, -0.14, 0.0},
    {0.0, -0.21, 0.0},
    {0.0, -0.28, 0.0},
    {0.0, -0.32, 0.0},   // sphere on link 4
    {0.0, -0.07, 0.0},
    {0.0, -0.14, 0.0},   // spheres on link 6         
    {-0.03, -0.06, -0.055},
    {0.0, -0.06, -0.055},
    {0.03, -0.06, -0.055}   // spheres on camera
};                                                            

// the radius of the spheres
const double SPHERE_RADIUS[NUM_SPHERES] = {0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.07, // spheres on link 2
                                           0.07, 0.06, 0.06, 0.06, 0.06,             // spheres on link 4
                                           0.05, 0.05,                               // spheres on link 6
                                           0.02, 0.02, 0.02};                        // spheres on camera   

// indices of the beginning and ending spheres of the tapered capsules for self-collision avoidance
// note that the last tapered capsule will be changed if gripper is considered (include_gripper_or_not = true)
const int NUM_TAPERED_CAPSULES = 3;
const std::pair<size_t, size_t> TAPERED_CAPSULES[NUM_TAPERED_CAPSULES] = {
    {0, 6},
    {7, 11},
    {12, 13}
}; 

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
                                const std::vector<Vec3>& boxCenters_input,
                                const std::vector<Vec3>& boxOrientation_input,
                                const std::vector<Vec3>& boxSize_input,
                                const bool include_gripper_or_not = false,
                                const double collision_buffer_input = 0,
                                Eigen::VectorXi jtype_input = Eigen::VectorXi(0));

    // Destructor
    ~KinovaCustomizedConstraints() = default;

    // class methods:
        // add one more sphere to the collision avoidance
    void add_collision_sphere(int joint_id,
                              const Vec3& offset,
                              double radius);

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

    std::unique_ptr<ForwardKinematicsSolver> fkPtr_;

    double collision_buffer = 0;

        // sphere info
    int num_spheres = 0;
    std::vector<int> sphere_joint_id;
    std::vector<Vec3> sphere_offset;
    std::vector<double> sphere_radius;

        // tapered capsule info
    std::vector<std::pair<size_t, size_t>> tapered_capsules;

    Eigen::Array<Vec3, Eigen::Dynamic, Eigen::Dynamic> sphere_centers_copy;
    Eigen::Array<MatX, Eigen::Dynamic, Eigen::Dynamic> sphere_centers_gradient_copy;

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
}; // namespace RAPTOR

#endif // KINOVA_CUSTOMIZED_CONSTRAINTS_H
