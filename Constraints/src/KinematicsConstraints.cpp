#include "KinematicsConstraints.h"

namespace IDTO {

KinematicsConstraints::KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                             const Model* model_input,
                                             const Eigen::VectorXi& jtype_input,
                                             const size_t joint_id_input,
                                             const size_t time_id_input,
                                             const Transform& desiredTransform_input,
                                             const Transform endT_input) :
    trajPtr_(trajPtr_input),
    modelPtr_(model_input),
    jtype(jtype_input),
    joint_id(joint_id_input),
    time_id(time_id_input),
    endT(endT_input) {
    fkPtr_ = std::make_unique<ForwardKinematicsSolver>(modelPtr_, jtype);

    if (joint_id > modelPtr_->nq) {
        throw std::invalid_argument("joint_id should not be larger than model.nq");
    }

    if (time_id >= trajPtr_->N) {
        throw std::invalid_argument("time_id should not be larger than number of instances in the trajectory");
    }

    desiredPosition = desiredTransform_input.p;
    desiredRotation = desiredTransform_input.R;

    if (Utils::ifTwoMatrixEqual(desiredRotation, -desiredRotation.transpose())) {
        throw std::invalid_argument("Input matrix is not skew-symmetric");
    }

    constrainPosition = true;
    constrainRotation = true;

    initialize_memory(6, trajPtr_->varLength);
}

KinematicsConstraints::KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                             const Model* model_input,
                                             const Eigen::VectorXi& jtype_input,
                                             const size_t joint_id_input,
                                             const size_t time_id_input,
                                             const Vec3& desiredPosition_input,
                                             const Transform endT_input) :
    trajPtr_(trajPtr_input),
    modelPtr_(model_input),
    jtype(jtype_input),
    joint_id(joint_id_input),
    time_id(time_id_input),
    endT(endT_input) {
    fkPtr_ = std::make_unique<ForwardKinematicsSolver>(modelPtr_, jtype);

    if (joint_id > modelPtr_->nq) {
        throw std::invalid_argument("joint_id should not be larger than model.nq");
    }

    if (time_id >= trajPtr_->N) {
        throw std::invalid_argument("time_id should not be larger than number of instances in the trajectory");
    }

    desiredPosition = desiredPosition_input;

    constrainPosition = true;
    constrainRotation = false;

    initialize_memory(3, trajPtr_->varLength);
}

KinematicsConstraints::KinematicsConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                             const Model* model_input,
                                             const Eigen::VectorXi& jtype_input,
                                             const size_t joint_id_input,
                                             const size_t time_id_input,
                                             const Mat3& desiredRotation_input,
                                             const Transform endT_input) :
    trajPtr_(trajPtr_input),
    modelPtr_(model_input),
    jtype(jtype_input),
    joint_id(joint_id_input),
    time_id(time_id_input),
    endT(endT_input) {
    fkPtr_ = std::make_unique<ForwardKinematicsSolver>(modelPtr_, jtype);

    if (joint_id > modelPtr_->nq) {
        throw std::invalid_argument("joint_id should not be larger than model.nq");
    }

    if (time_id >= trajPtr_->N) {
        throw std::invalid_argument("time_id should not be larger than number of instances in the trajectory");
    }

    desiredRotation = desiredRotation_input;
    
    if (Utils::ifTwoMatrixEqual(desiredRotation, -desiredRotation.transpose())) {
        throw std::invalid_argument("Input matrix is not skew-symmetric");
    }

    constrainPosition = false;
    constrainRotation = true;

    initialize_memory(3, trajPtr_->varLength);
}

void KinematicsConstraints::compute(const VecX& z, 
                                    bool compute_derivatives,
                                    bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);
    
    const VecX& q = trajPtr_->q(time_id);
    const MatX& pq_pz = trajPtr_->pq_pz(time_id);
    const Eigen::Array<MatX, Eigen::Dynamic, 1>& pq_pz_pz = trajPtr_->pq_pz_pz.col(time_id);

    if (compute_hessian) {
        fkPtr_->compute(0, joint_id, q, &startT, &endT, 2);
    }
    else if (compute_derivatives) {
        fkPtr_->compute(0, joint_id, q, &startT, &endT, 1);
    }
    else {
        fkPtr_->compute(0, joint_id, q, &startT, &endT, 0);
    }

    MatX pg_pq;
    Eigen::Array<MatX, 3, 1> pg_pq_pq;

    if (constrainPosition) {
        if (compute_hessian) {
            g.head(3) = LieSpaceResidual::translationResidual(fkPtr_, desiredPosition, &pg_pq, &pg_pq_pq);
            pg_pz.topRows(3) = pg_pq * pq_pz;
            for (int i = 0; i < 3; i++) {
                pg_pz_pz(i) = pq_pz.transpose() * pg_pq_pq(i) * pq_pz;

                for (int j = 0; j < pq_pz_pz.size(); j++) {
                    pg_pz_pz(i) += pg_pq(i, j) * pq_pz_pz(j);
                }
            }

        }
        else if (compute_derivatives) {
            g.head(3) = LieSpaceResidual::translationResidual(fkPtr_, desiredPosition, &pg_pq);
            pg_pz.topRows(3) = pg_pq * pq_pz;
        }
        else {
            g.head(3) = LieSpaceResidual::translationResidual(fkPtr_, desiredPosition);
        }
    }
    
    if (constrainRotation) {
        if (compute_hessian) {
            g.tail(3) = LieSpaceResidual::rotationResidual(fkPtr_, desiredRotation, &pg_pq, &pg_pq_pq);
            pg_pz.bottomRows(3) = pg_pq * pq_pz;

            int startIdx = constrainPosition ? 3 : 0;

            for (int i = startIdx; i < startIdx + 3; i++) {
                pg_pz_pz(i) = pq_pz.transpose() * pg_pq_pq(i - startIdx) * pq_pz;

                for (int j = 0; j < pq_pz_pz.size(); j++) {
                    pg_pz_pz(i) += pg_pq(i - startIdx, j) * pq_pz_pz(j);
                }
            }
        }
        else if (compute_derivatives) {
            g.tail(3) = LieSpaceResidual::rotationResidual(fkPtr_, desiredRotation, &pg_pq);
            pg_pz.bottomRows(3) = pg_pq * pq_pz;
        }
        else {
            g.tail(3) = LieSpaceResidual::rotationResidual(fkPtr_, desiredRotation);
        }
    }
}

void KinematicsConstraints::compute_bounds() {
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
}

void KinematicsConstraints::print_violation_info() {
    if (constrainPosition) {
        const VecX& g_position = g.head(3);

        if (abs(g_position(0)) > 1e-5) {
            std::cout << "    Error on position x: " 
                      << g_position(0) 
                      << std::endl;
        }
        if (abs(g_position(1)) > 1e-5) {
            std::cout << "    Error on position y: " 
                      << g_position(1) << std::endl;
        }
        if (abs(g_position(2)) > 1e-5) {
            std::cout << "    Error on position z: " 
                      << g_position(2) << std::endl;
        }
    }

    if (constrainRotation) {
        const VecX& g_rotation = g.tail(3);

        if (g_rotation.norm() > 1e-5) {
            std::cout << "    Error on rotation (norm of skew residual): " 
                      << g_rotation.norm() << std::endl;
        }
    }
}

}; // namespace IDTO