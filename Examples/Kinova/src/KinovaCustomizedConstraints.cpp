#include "KinovaCustomizedConstraints.h"

namespace IDTO {
namespace Kinova {

KinovaCustomizedConstraints::KinovaCustomizedConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                                         const Model& model_input,
                                                         const Eigen::VectorXi& jtype_input,
                                                         const Eigen::Array<Vec3, 1, Eigen::Dynamic>& boxCenters_input,
                                                         const Eigen::Array<Vec3, 1, Eigen::Dynamic>& boxOrientation_input,
                                                         const Eigen::Array<Vec3, 1, Eigen::Dynamic>& boxSize_input) :
    trajPtr_(trajPtr_input),
    jtype(jtype_input) {
    modelPtr_ = std::make_unique<Model>(model_input);
    fkhofPtr_ = std::make_unique<ForwardKinematicsHighOrderDerivative>();
    collisionAvoidancePtr_ = std::make_shared<BoxCollisionAvoidance>(boxCenters_input, 
                                                                     boxOrientation_input,
                                                                     boxSize_input);

    // m = trajPtr_->N * NUM_SPHERES * collisionAvoidancePtr_->numObstacles;
    m = trajPtr_->N * NUM_SPHERES;

    jointTJ = MatX::Zero(3, trajPtr_->Nact);

    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, trajPtr_->varLength);
}

void KinovaCustomizedConstraints::compute(const VecX& z, 
                                          bool compute_derivatives,
                                          bool compute_hessian) {
    if (compute_hessian) {
        throw std::invalid_argument("KinovaCustomizedConstraints does not support hessian computation");
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);
    
    const int numObstacles = collisionAvoidancePtr_->numObstacles;

    if (numObstacles == 0) {
        return;
    }

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

            // g.block((i * NUM_SPHERES + j) * numObstacles, 
            //         0, 
            //         numObstacles, 
            //         1) = collisionAvoidancePtr_->distances.array() - sphere_radius[j];

            // if (compute_derivatives) {
            //     pg_pz.block((i * NUM_SPHERES + j) * numObstacles, 
            //                 0, 
            //                 numObstacles, 
            //                 trajPtr_->varLength) = collisionAvoidancePtr_->pdistances_pz;
            // }
            Eigen::VectorXd::Index minIndex;
            collisionAvoidancePtr_->distances.minCoeff(&minIndex);
            double distance = collisionAvoidancePtr_->distances(minIndex) - sphere_radius[j];
            g(i * NUM_SPHERES + j) = distance;

            if (compute_derivatives) {
                pg_pz.row(i * NUM_SPHERES + j) = collisionAvoidancePtr_->pdistances_pz.row(minIndex);
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
    const int numObstacles = collisionAvoidancePtr_->numObstacles;
    // for (int i = 0; i < trajPtr_->N; i++) {
    //     for (int j = 0; j < NUM_SPHERES; j++) {
    //         for (int k = 0; k < numObstacles; k++) {
    //             if (g((i * NUM_SPHERES + j) * numObstacles + k) <= 0) {
    //                 std::cout << "        KinovaCustomizedConstraints.cpp: Sphere " << j 
    //                           << " corresponding to link " << sphere_joint_id[j]
    //                           << " at time instance " << i 
    //                           << " is in collision with obstacle " << k
    //                           << std::endl;
    //             }
    //         }
    //     }
    // }
    for (int i = 0; i < trajPtr_->N; i++) {
        for (int j = 0; j < NUM_SPHERES; j++) {
            if (g(i * NUM_SPHERES + j) <= 0) {
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
}; // namespace IDTO