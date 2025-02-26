#include "ConstrainedInverseDynamics.h"

namespace RAPTOR {

ConstrainedInverseDynamics::ConstrainedInverseDynamics(const Model& model_input, 
                                                       std::shared_ptr<Trajectories>& trajPtr_input,
                                                       int numDependentJoints_input) : 
    InverseDynamics(model_input, trajPtr_input),
    numDependentJoints(numDependentJoints_input) {
    numIndependentJoints = modelPtr_->nv - numDependentJoints;

    tau_dep = VecX::Zero(numDependentJoints);
    tau_indep = VecX::Zero(numIndependentJoints);

    lambda.resize(1, N);
    plambda_pz.resize(1, N);

    q.resize(1, N);
    v.resize(1, N);
    a.resize(1, N);

    pq_pz.resize(1, N);
    pv_pz.resize(1, N);
    pa_pz.resize(1, N);

    prnea_pq_dep = MatX::Zero(numDependentJoints, modelPtr_->nv);
    prnea_pq_indep = MatX::Zero(numIndependentJoints, modelPtr_->nv);
    prnea_pv_dep = MatX::Zero(numDependentJoints, modelPtr_->nv);
    prnea_pv_indep = MatX::Zero(numIndependentJoints, modelPtr_->nv);
    prnea_pa_dep = MatX::Zero(numDependentJoints, modelPtr_->nv);
    prnea_pa_indep = MatX::Zero(numIndependentJoints, modelPtr_->nv);

    plambda_pq = MatX::Zero(numDependentJoints, modelPtr_->nv);
    plambda_pv = MatX::Zero(numDependentJoints, modelPtr_->nv);
    plambda_pa = MatX::Zero(numDependentJoints, modelPtr_->nv);

    JTlambda_partial_dq = MatX::Zero(modelPtr_->nv, modelPtr_->nv);
    JTlambda_partial_dq_dep = MatX::Zero(numDependentJoints, modelPtr_->nv);
    JTlambda_partial_dq_indep = MatX::Zero(numIndependentJoints, modelPtr_->nv);

    // dcPtr_ is still empty at this moment!!!
    // You need to define it in the constructor of your own derived class!!!
}

void ConstrainedInverseDynamics::compute(const VecX& z,
                                         bool compute_derivatives,
                                         bool compute_hessian) {
    if (trajPtr_ == nullptr) {
        throw std::runtime_error("trajPtr_ is not defined yet!");
    }

    if (dcPtr_ == nullptr) {
        throw std::runtime_error("dcPtr_ is not defined yet!");
    }                 

    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }         

    if (compute_hessian) {
        throw std::invalid_argument("ConstrainedInverseDynamics does not support hessian computation");
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    for (int i = 0; i < N; i++) {  
        // fill in independent indeces in a vector
        q(i) = VecX::Zero(modelPtr_->nv);
        v(i) = VecX::Zero(modelPtr_->nv);
        a(i) = VecX::Zero(modelPtr_->nv);
        dcPtr_->fill_independent_vector(q(i), trajPtr_->q(i));
        dcPtr_->fill_independent_vector(v(i), trajPtr_->q_d(i));
        dcPtr_->fill_independent_vector(a(i), trajPtr_->q_dd(i));
                           
        // always call this first to update dcPtr_->J_dep and dcPtr_->J_indep!!!   
        dcPtr_->setupJointPositionVelocityAcceleration(q(i), v(i), a(i), compute_derivatives);

        if (compute_derivatives) {
            pq_pz(i) = MatX::Zero(modelPtr_->nv, trajPtr_->varLength);
            pv_pz(i) = MatX::Zero(modelPtr_->nv, trajPtr_->varLength);
            pa_pz(i) = MatX::Zero(modelPtr_->nv, trajPtr_->varLength);

            // fill in independent joints derivatives directly
            for (int j = 0; j < dcPtr_->numIndependentJoints; j++) {
                int indenpendentJointIndex = dcPtr_->return_independent_joint_index(j);
                pq_pz(i).row(indenpendentJointIndex) = trajPtr_->pq_pz(i).row(j);
                pv_pz(i).row(indenpendentJointIndex) = trajPtr_->pq_d_pz(i).row(j);
                pa_pz(i).row(indenpendentJointIndex) = trajPtr_->pq_dd_pz(i).row(j);
            }

            // derivatives of unactuated joint positions
            MatX pq_dep_pz = dcPtr_->pq_dep_pq_indep * trajPtr_->pq_pz(i);

            for (int j = 0; j < dcPtr_->numDependentJoints; j++) {
                int denpendentJointIndex = dcPtr_->return_dependent_joint_index(j);
                pq_pz(i).row(denpendentJointIndex) = pq_dep_pz.row(j);
            }

            // derivatives of unactuated joint velocities
            MatX pv_dep_pz = dcPtr_->pv_dep_pq       * pq_pz(i) +            // all joint positions  
                             dcPtr_->pv_dep_pv_indep * trajPtr_->pq_d_pz(i); // actuated joint velocities

            for (int j = 0; j < dcPtr_->numDependentJoints; j++) {
                int denpendentJointIndex = dcPtr_->return_dependent_joint_index(j);
                pv_pz(i).row(denpendentJointIndex) = pv_dep_pz.row(j);
            }

            // derivatives of unactuated joint accelerations
            MatX pa_dep_pz = dcPtr_->pa_dep_pq       * pq_pz(i) +             // all joint positions  
                             dcPtr_->pa_dep_pv       * pv_pz(i) +             // all joint velocities  
                             dcPtr_->pa_dep_pa_indep * trajPtr_->pq_dd_pz(i); // actuated joint accelerations

            for (int j = 0; j < dcPtr_->numDependentJoints; j++) {
                int denpendentJointIndex = dcPtr_->return_dependent_joint_index(j);
                pa_pz(i).row(denpendentJointIndex) = pa_dep_pz.row(j);
            }
        }

        // compute unconstrained inverse dynamics
        if (!compute_derivatives) {
            pinocchio::rnea(*modelPtr_, *dataPtr_, 
                            q(i), v(i), a(i));
        }
        else {
            pinocchio::computeRNEADerivatives(*modelPtr_, *dataPtr_, 
                                              q(i), v(i), a(i), 
                                              prnea_pq, prnea_pv, prnea_pa);
        }

        // adjust with damping force
        tau(i) = dataPtr_->tau + 
                 modelPtr_->damping.cwiseProduct(v(i));
                
        if (compute_derivatives) {
            // prnea_pa is just the inertia matrix.
            // pinocchio only computes the upper triangle part of it.
            for (int mi = 0; mi < prnea_pa.rows(); mi++) {
                for (int mj = 0; mj <= mi; mj++) {
                    prnea_pa(mi, mj) = prnea_pa(mj, mi);
                }
            }

            // pinocchio rnea does not take damping into account
            for (int mi = 0; mi < prnea_pa.rows(); mi++) {
                prnea_pv(mi, mi) += modelPtr_->damping(mi);
            }
        }

        tau_dep = dcPtr_->get_dependent_vector(tau(i));
        tau_indep = dcPtr_->get_independent_vector(tau(i));

        // assume setupJointPositionVelocityAcceleration() has been called
        // compute the reaction force to maintain constraints
        lambda(i) = dcPtr_->J_dep_T_qr.solve(tau_dep);

        if (compute_derivatives) {
            dcPtr_->get_dependent_rows(prnea_pq_dep, prnea_pq);
            dcPtr_->get_independent_rows(prnea_pq_indep, prnea_pq);
            dcPtr_->get_dependent_rows(prnea_pv_dep, prnea_pv);
            dcPtr_->get_independent_rows(prnea_pv_indep, prnea_pv);
            dcPtr_->get_dependent_rows(prnea_pa_dep, prnea_pa);
            dcPtr_->get_independent_rows(prnea_pa_indep, prnea_pa);

            dcPtr_->get_JTx_partial_dq(q(i), lambda(i));
            JTlambda_partial_dq = dcPtr_->JTx_partial_dq;
            dcPtr_->get_dependent_rows(JTlambda_partial_dq_dep, JTlambda_partial_dq);
            dcPtr_->get_independent_rows(JTlambda_partial_dq_indep, JTlambda_partial_dq);

            plambda_pz(i) = dcPtr_->J_dep_T_qr.solve((prnea_pq_dep - JTlambda_partial_dq_dep) * pq_pz(i) + 
                                                      prnea_pv_dep * pv_pz(i) + 
                                                      prnea_pa_dep * pa_pz(i));
        }

        // assume setupJointPositionVelocityAcceleration() has been called
        tau(i) = tau_indep - dcPtr_->J_indep.transpose() * lambda(i);

        if (compute_derivatives) {
            ptau_pz(i) = (prnea_pq_indep - JTlambda_partial_dq_indep) * pq_pz(i) + 
                         prnea_pv_indep * pv_pz(i) + 
                         prnea_pa_indep * pa_pz(i) - 
                         dcPtr_->J_indep.transpose() * plambda_pz(i);
        }
    }
}

}; // namespace RAPTOR

