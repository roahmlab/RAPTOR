#include "ConstrainedInverseDynamics.h"

namespace IDTO {

ConstrainedInverseDynamics::ConstrainedInverseDynamics(const Model& model_input, 
                                                       int N_input, 
                                                       int numDependentJoints_input) : 
    InverseDynamics(model_input, N_input),
    numDependentJoints(numDependentJoints_input) {
    numIndependentJoints = modelPtr_->nv - numDependentJoints;

    tau_dep = VecX::Zero(numDependentJoints);
    tau_indep = VecX::Zero(numIndependentJoints);

    lambda.resize(1, N);

    plambda_pq.resize(1, N);
    plambda_pv.resize(1, N);
    plambda_pa.resize(1, N);

    rnea_partial_dq_dep = MatX::Zero(numDependentJoints, modelPtr_->nv);
    rnea_partial_dq_indep = MatX::Zero(numIndependentJoints, modelPtr_->nv);
    rnea_partial_dv_dep = MatX::Zero(numDependentJoints, modelPtr_->nv);
    rnea_partial_dv_indep = MatX::Zero(numIndependentJoints, modelPtr_->nv);
    rnea_partial_da_dep = MatX::Zero(numDependentJoints, modelPtr_->nv);
    rnea_partial_da_indep = MatX::Zero(numIndependentJoints, modelPtr_->nv);

    JTlambda_partial_dq = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    JTlambda_partial_dq_dep = MatX::Zero(numDependentJoints, modelPtr_->nv);
    JTlambda_partial_dq_indep = MatX::Zero(numIndependentJoints, modelPtr_->nv);

    // dynamicsConstraintsPtr_ is still empty at this moment!!!
    // You need to define it in the constructor of your own derived class!!!
}

void ConstrainedInverseDynamics::compute(Eigen::Array<VecX, 1, Eigen::Dynamic>& q, 
                                         Eigen::Array<VecX, 1, Eigen::Dynamic>& v, 
                                         Eigen::Array<VecX, 1, Eigen::Dynamic>& a,
                                         bool compute_derivatives) {
    for (int i = 0; i < N; i++) {                     
        // always call this first to update dynamicsConstraintsPtr_->J_dep and dynamicsConstraintsPtr_->J_indep!!!   
        dynamicsConstraintsPtr_->setupJointPositionVelocityAcceleration(q(i), v(i), a(i), compute_derivatives);

        if (!compute_derivatives) {
            pinocchio::rnea(*modelPtr_, *dataPtr_, q(i), v(i), a(i));
        }
        else {
            pinocchio::computeRNEADerivatives(*modelPtr_, *dataPtr_, q(i), v(i), a(i), 
                                              rnea_partial_dq, rnea_partial_dv, rnea_partial_da);
        }

        // adjust with damping force and rotor inertia force
        tau(i) = dataPtr_->tau + 
                 modelPtr_->damping.cwiseProduct(v(i)) + 
                 modelPtr_->rotorInertia.cwiseProduct(a(i));
                
        if (compute_derivatives) {
            // rnea_partial_da is just the inertia matrix.
            // pinocchio only computes the upper triangle part of it.
            for (int mi = 0; mi < rnea_partial_da.rows(); mi++) {
                for (int mj = 0; mj <= mi; mj++) {
                    if (mi == mj) {
                        rnea_partial_da(mi, mj) += modelPtr_->rotorInertia(mi);
                    }
                    else {
                        rnea_partial_da(mi, mj) = rnea_partial_da(mj, mi);
                    }
                }
            }

            // pinocchio rnea does not take damping into account
            for (int mi = 0; mi < rnea_partial_da.rows(); mi++) {
                rnea_partial_dv(mi, mi) += modelPtr_->damping(mi);
            }
        }

        tau_dep = dynamicsConstraintsPtr_->get_dependent_vector(tau(i));
        tau_indep = dynamicsConstraintsPtr_->get_independent_vector(tau(i));

        // assume setupJointPositionVelocityAcceleration() has been called
        lambda(i) = dynamicsConstraintsPtr_->J_dep_T_qr.solve(tau_dep);

        if (compute_derivatives) {
            dynamicsConstraintsPtr_->get_dependent_rows(rnea_partial_dq_dep, rnea_partial_dq);
            dynamicsConstraintsPtr_->get_independent_rows(rnea_partial_dq_indep, rnea_partial_dq);
            dynamicsConstraintsPtr_->get_dependent_rows(rnea_partial_dv_dep, rnea_partial_dv);
            dynamicsConstraintsPtr_->get_independent_rows(rnea_partial_dv_indep, rnea_partial_dv);
            dynamicsConstraintsPtr_->get_dependent_rows(rnea_partial_da_dep, rnea_partial_da);
            dynamicsConstraintsPtr_->get_independent_rows(rnea_partial_da_indep, rnea_partial_da);

            dynamicsConstraintsPtr_->get_JTx_partial_dq(q(i), lambda(i));
            MatX JTlambda_partial_dq = dynamicsConstraintsPtr_->JTx_partial_dq;
            dynamicsConstraintsPtr_->get_dependent_rows(JTlambda_partial_dq_dep, JTlambda_partial_dq);
            dynamicsConstraintsPtr_->get_independent_rows(JTlambda_partial_dq_indep, JTlambda_partial_dq);

            plambda_pq(i) = dynamicsConstraintsPtr_->J_dep_T_qr.solve(rnea_partial_dq_dep - JTlambda_partial_dq_dep);
            plambda_pv(i) = dynamicsConstraintsPtr_->J_dep_T_qr.solve(rnea_partial_dv_dep);
            plambda_pa(i) = dynamicsConstraintsPtr_->J_dep_T_qr.solve(rnea_partial_da_dep);
        }

        // assume setupJointPositionVelocityAcceleration() has been called
        tau(i) = tau_indep - dynamicsConstraintsPtr_->J_indep.transpose() * lambda(i);

        if (compute_derivatives) {
            ptau_pq(i) = rnea_partial_dq_indep - dynamicsConstraintsPtr_->J_indep.transpose() * plambda_pq(i) - JTlambda_partial_dq_indep;
            ptau_pv(i) = rnea_partial_dv_indep - dynamicsConstraintsPtr_->J_indep.transpose() * plambda_pv(i);
            ptau_pa(i) = rnea_partial_da_indep - dynamicsConstraintsPtr_->J_indep.transpose() * plambda_pa(i);
        }
    }
}

}; // namespace IDTO

