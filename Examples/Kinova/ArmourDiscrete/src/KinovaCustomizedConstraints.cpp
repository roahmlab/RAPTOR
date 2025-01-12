#include "KinovaCustomizedConstraints.h"

namespace RAPTOR {
namespace Kinova {

KinovaCustomizedConstraints::KinovaCustomizedConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                                         const Model& model_input,
                                                         const std::vector<Vec3>& boxCenters_input,
                                                         const std::vector<Vec3>& boxOrientation_input,
                                                         const std::vector<Vec3>& boxSize_input,
                                                         const bool include_gripper_or_not,
                                                         const double collision_buffer_input,
                                                         Eigen::VectorXi jtype_input) :
    trajPtr_(trajPtr_input),
    collision_buffer(collision_buffer_input) {
    modelPtr_ = std::make_unique<Model>(model_input);
    fkPtr_ = std::make_unique<ForwardKinematicsSolver>(modelPtr_.get(), jtype_input);
    collisionAvoidancePtr_ = std::make_shared<BoxCollisionAvoidance>(boxCenters_input, 
                                                                     boxOrientation_input,
                                                                     boxSize_input);

    // initialize sphere info
    num_spheres = NUM_SPHERES;
    sphere_joint_id.resize(NUM_SPHERES);
    sphere_offset.resize(NUM_SPHERES);
    sphere_radius.resize(NUM_SPHERES);
    for (int i = 0; i < NUM_SPHERES; i++) {
        sphere_joint_id[i] = SPHERE_JOINT_ID[i];
        sphere_offset[i] = Vec3(SPHERE_OFFSET[i][0], 
                                SPHERE_OFFSET[i][1], 
                                SPHERE_OFFSET[i][2]);
        sphere_radius[i] = SPHERE_RADIUS[i];
    }

    // initialize tapered capsules
    tapered_capsules.reserve(NUM_TAPERED_CAPSULES);
    for (int i = 0; i < NUM_TAPERED_CAPSULES; i++) {
        tapered_capsules.push_back(TAPERED_CAPSULES[i]);
    }

    // 3 spheres for gripper
    if (include_gripper_or_not) {
        add_collision_sphere(7, Vec3(0.0, 0.0, -0.10), 0.05); // sphere # 17
        add_collision_sphere(7, Vec3(0.0, 0.0, -0.15), 0.05); // sphere # 18
        add_collision_sphere(7, Vec3(0.0, 0.0, -0.20), 0.05); // sphere # 19
        tapered_capsules.back().second = 19; // the last tapered capsule is now extended to the gripper
    }

    // m = trajPtr_->N * NUM_SPHERES * collisionAvoidancePtr_->numObstacles;

    collisionAvoidancePtr_->onlyComputeDerivativesForMinimumDistance = true;
    m = trajPtr_->N * num_spheres;

    jointTJ = MatX::Zero(3, trajPtr_->Nact);
    sphere_centers_copy.resize(num_spheres, trajPtr_->N);
    sphere_centers_gradient_copy.resize(num_spheres, trajPtr_->N);

    initialize_memory(m, trajPtr_->varLength);
}

void KinovaCustomizedConstraints::add_collision_sphere(int joint_id,
                                                       const Vec3& offset,
                                                       double radius) {
    num_spheres++;
    sphere_joint_id.push_back(joint_id);
    sphere_offset.push_back(offset);
    sphere_radius.push_back(radius);

    m = trajPtr_->N * num_spheres;
    initialize_memory(m, trajPtr_->varLength);
}

void KinovaCustomizedConstraints::compute(const VecX& z, 
                                          bool compute_derivatives,
                                          bool compute_hessian) {
    trajPtr_->compute(z, compute_derivatives, compute_hessian);
    
    const int numObstacles = collisionAvoidancePtr_->numObstacles;

    if (numObstacles == 0) {
        return;
    }

    Vec3 sphereCenter;
    MatX psphereCenter_pz;
    Eigen::Array<MatX, 3, 1> psphereCenter_pz_pz;

    for (int i = 0; i < trajPtr_->N; i++) {
        const VecX& q = trajPtr_->q(i).head(trajPtr_->Nact);
        const MatX& pq_pz = trajPtr_->pq_pz(i);
        const Eigen::Array<MatX, Eigen::Dynamic, 1>& pq_pz_pz = trajPtr_->pq_pz_pz.col(i); 

        for (int j = 0; j < num_spheres; j++) {
            // define the transform matrix of the sphere center with respect to the joint
            endT = Transform(sphere_offset[j]);

            if (compute_hessian) {
                fkPtr_->compute(0, sphere_joint_id[j], q, &startT, &endT, 2);
            }
            else if (compute_derivatives) {
                fkPtr_->compute(0, sphere_joint_id[j], q, &startT, &endT, 1);
            }
            else {
                fkPtr_->compute(0, sphere_joint_id[j], q, &startT, &endT, 0);
            }

            sphereCenter = fkPtr_->getTranslation();

            sphere_centers_copy(j, i) = sphereCenter;

            if (compute_derivatives) {
                psphereCenter_pz = fkPtr_->getTranslationJacobian() * pq_pz;

                sphere_centers_gradient_copy(j, i) = psphereCenter_pz;
            }

            if (compute_hessian) {
                Eigen::Array<MatX, 3, 1> psphereCenter_pq_pq;
                fkPtr_->getTranslationHessian(psphereCenter_pq_pq);

                // (1) p2_FK_pq2 * pq_pz1 * pq_pz2
                for (int k = 0; k < 3; k++) {
                    psphereCenter_pz_pz(k) = pq_pz.transpose() * psphereCenter_pq_pq(k) * pq_pz;
                }
                
                // (2) p_FK_pq * p2_q_pz2 (this should be zero for most trajectories)
                for (int k = 0; k < 3; k++) {
                    for (int p = 0; p < trajPtr_->Nact; p++) {
                        psphereCenter_pz_pz(k) += jointTJ(k, p) * pq_pz_pz(p);
                    }
                }
            }

            if (compute_hessian) {
                collisionAvoidancePtr_->computeDistance(sphereCenter, psphereCenter_pz, psphereCenter_pz_pz);
            }
            else if (compute_derivatives) {
                collisionAvoidancePtr_->computeDistance(sphereCenter, psphereCenter_pz);
            }
            else {
                collisionAvoidancePtr_->computeDistance(sphereCenter);
            }
            
            const auto& minIndex = collisionAvoidancePtr_->minimumDistanceIndex;
            g(i * num_spheres + j) = collisionAvoidancePtr_->minimumDistance - sphere_radius[j];

            if (compute_derivatives) {
                pg_pz.row(i * num_spheres + j) = 
                    collisionAvoidancePtr_->pdistances_pz.row(minIndex);
            }

            if (compute_hessian) {
                pg_pz_pz(i * num_spheres + j) = 
                    collisionAvoidancePtr_->pdistances_pz_pz(minIndex);
            }
        }
    }
}

void KinovaCustomizedConstraints::compute_bounds() {
    // distance with obstacles larger than collision_buffer = 0
    g_lb.setConstant(collision_buffer);
    g_ub.setConstant(1e19);
}

void KinovaCustomizedConstraints::print_violation_info() {
    for (int i = 0; i < trajPtr_->N; i++) {
        for (int j = 0; j < num_spheres; j++) {
            if (g(i * num_spheres + j) <= 0) {
                std::cout << "        KinovaCustomizedConstraints.cpp: Sphere " << j 
                          << " corresponding to link " << sphere_joint_id[j]
                          << " at time instance " << i 
                          << " is in collision with the environment"
                          << std::endl;
            }
        }
    }
}

}; // namespace Kinova
}; // namespace RAPTOR