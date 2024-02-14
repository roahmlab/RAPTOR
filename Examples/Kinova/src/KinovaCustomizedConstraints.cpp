#include "KinovaCustomizedConstraints.h"

namespace IDTO {
namespace Kinova {

KinovaCustomizedConstraints::KinovaCustomizedConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                                         const Model& model_input,
                                                         const Eigen::VectorXi& jtype_input,
                                                         const Eigen::Array<Vec3, 1, Eigen::Dynamic>& zonotopeCenters_input,
                                                         const Eigen::Array<MatX, 1, Eigen::Dynamic>& zonotopeGenerators_input) :
    trajPtr_(trajPtr_input),
    jtype(jtype_input) {
    modelPtr_ = std::make_unique<Model>(model_input);
    fkhofPtr_ = std::make_unique<ForwardKinematicsHighOrderDerivative>();
    collisionAvoidancePtr_ = std::make_shared<ZonotopeCollisionAvoidance>(zonotopeCenters_input, 
                                                                          zonotopeGenerators_input);

    m = trajPtr_->N * NUM_SPHERES * collisionAvoidancePtr_->numObstacles;

    jointTJ = MatX::Zero(3, trajPtr_->Nact);

    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, trajPtr_->varLength);
}

void KinovaCustomizedConstraints::compute(const VecX& z, bool compute_derivatives) {
    trajPtr_->compute(z, compute_derivatives);
    
    for (int i = 0; i < trajPtr_->N; i++) {
        const VecX& q = trajPtr_->q(i).head(trajPtr_->Nact);

        for (int j = 0; j < NUM_SPHERES; j++) {
            // define the transform matrix of the sphere center with respect to the joint
            endT = Transform(sphere_offset_axis[j], sphere_offset[j]);

            fkhofPtr_->fk(jointT, *modelPtr_, jtype, sphere_joint_id[j], 0, q, endT, startT);

            Vec3 sphereCenters = fkhofPtr_->Transform2xyz(jointT);

            collisionAvoidancePtr_->computeDistance(sphereCenters);

            if (compute_derivatives) {
                fkhofPtr_->fk_jacobian(dTdq, *modelPtr_, jtype, sphere_joint_id[j], 0, q, endT, startT);
                fkhofPtr_->Transform2xyzJacobian(jointTJ, jointT, dTdq);

                MatX psphereCenters_pz = jointTJ * trajPtr_->pq_pz(i);

                collisionAvoidancePtr_->computeDistance(sphereCenters, psphereCenters_pz);
            }

            g.block((i * NUM_SPHERES + j) * collisionAvoidancePtr_->numObstacles, 
                    0, 
                    collisionAvoidancePtr_->numObstacles, 
                    1) = collisionAvoidancePtr_->distances.array() - sphere_radius[j];

            if (compute_derivatives) {
                pg_pz.block((i * NUM_SPHERES + j) * collisionAvoidancePtr_->numObstacles, 
                            0, 
                            collisionAvoidancePtr_->numObstacles, 
                            trajPtr_->varLength) = collisionAvoidancePtr_->pdistances_pz;
            }
        }
    }
}

void KinovaCustomizedConstraints::compute_bounds() {
    // distance with obstacles larger than 0
    g_lb.setConstant(0);
    g_ub.setConstant(1e19);
}

void KinovaCustomizedConstraints::print_violation_info() {
    for (int i = 0; i < trajPtr_->N; i++) {
        for (int j = 0; j < NUM_SPHERES; j++) {
            for (int k = 0; k < collisionAvoidancePtr_->numObstacles; k++) {
                if (g((i * NUM_SPHERES + j) * collisionAvoidancePtr_->numObstacles + k) <= 0) {
                    std::cout << "        KinovaCustomizedConstraints.cpp: Sphere " 
                              << j 
                              << " corresponding to link "
                              << sphere_joint_id[j]
                              << " at time instance " 
                              << i 
                              << " is in collision with obstacle " 
                              << k
                              << std::endl;
                }
            }
        }
    }
}

}; // namespace Kinova
}; // namespace IDTO