#include "KinovaCustomizedConstraints.h"

namespace RAPTOR {
namespace Kinova {

KinovaCustomizedConstraints::KinovaCustomizedConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                                         const Model& model_input,
                                                         const Eigen::VectorXi& jtype_input,
                                                         const std::vector<Vec3>& boxCenters_input,
                                                         const std::vector<Vec3>& boxOrientation_input,
                                                         const std::vector<Vec3>& boxSize_input) :
    trajPtr_(trajPtr_input),
    jtype(jtype_input) {
    modelPtr_ = std::make_unique<Model>(model_input);
    fkPtr_ = std::make_unique<ForwardKinematicsSolver>(modelPtr_.get(), jtype);
    collisionAvoidancePtr_ = std::make_shared<BoxCollisionAvoidance>(boxCenters_input, 
                                                                     boxOrientation_input,
                                                                     boxSize_input);

    // m = trajPtr_->N * NUM_SPHERES * collisionAvoidancePtr_->numObstacles;

    collisionAvoidancePtr_->onlyComputeDerivativesForMinimumDistance = true;
    m = trajPtr_->N * NUM_SPHERES;

    jointTJ = MatX::Zero(3, trajPtr_->Nact);

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

        for (int j = 0; j < NUM_SPHERES; j++) {
            // define the transform matrix of the sphere center with respect to the joint
            endT = Transform(sphere_offset_axis[j], sphere_offset[j]);

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

            if (compute_derivatives) {
                psphereCenter_pz = fkPtr_->getTranslationJacobian() * pq_pz;
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
            g(i * NUM_SPHERES + j) = collisionAvoidancePtr_->minimumDistance - sphere_radius[j];

            if (compute_derivatives) {
                pg_pz.row(i * NUM_SPHERES + j) = 
                    collisionAvoidancePtr_->pdistances_pz.row(minIndex);
            }

            if (compute_hessian) {
                pg_pz_pz(i * NUM_SPHERES + j) = 
                    collisionAvoidancePtr_->pdistances_pz_pz(minIndex);
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
}; // namespace RAPTOR