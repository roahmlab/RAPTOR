#include "WaitrContactConstraints.h"

namespace IDTO {

WaitrContactConstraints::WaitrContactConstraints(std::shared_ptr<CustomizedInverseDynamics>& idPtr_input,
                                                 const contactSurfaceParams& csp_input) :
    idPtr_(idPtr_input),
    fp(csp_input) {
    // (1)    positive contact force
    // (2)    translation friction cone
    // (3)    rotation friction cone
    // (4, 5) ZMP on one axis
    // (6, 7) ZMP on the other axis
    m = idPtr_->trajPtr_->N * 6; // TODO: this is the number of constraints for all time instances

    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, idPtr_->trajPtr_->varLength);
}

void WaitrContactConstraints::compute(const VecX& z, bool compute_derivatives) {
    if (is_computed(z, compute_derivatives)) {
        return;
    }

    if (compute_derivatives) {
        pg_pz.setZero();
    }
    
    idPtr_->compute(z, compute_derivatives);

    for (int i = 0; i < idPtr_->trajPtr_->N; i++) {
        // access the last contact wrench, 
        // which should be the contact wrench between the end effector and the object
        const Vec6& lambda = idPtr_->lambda(i);

        const Vec3& rotation_torque = lambda.head(3); // This is a 3D vector
        const Vec3& translation_force = lambda.tail(3); // This is a 3D vector

        // TODO: edit the contact constraints here for each time instance, 
        // based on translation_force and rotation_torque
        // You can refer to SurfaceContactConstraints.cpp, which should be very similar
        g.block(6 * i, 0, 3, 1) = rotation_torque;
        g.block(6 * i + 3, 0, 3, 1) = translation_force;

        if (compute_derivatives) {
            const MatX& plambda_pz = idPtr_->plambda_pz(i);

            // TODO: edit the gradient
            // plambda_pz.topRows(3) will be the gradient of the translation force w.r.t z
            // plambda_pz.bottomRows(3) will be the gradient of the rotation torque w.r.t z
            pg_pz.block(6 * i, 0, 6, idPtr_->trajPtr_->varLength) = plambda_pz;
        }
    }
}

void WaitrContactConstraints::compute_bounds() {
    // TODO: fill in the lower and upper bounds for the constraints
    // This is consistent with what ipopt is looking for
    // g_lb and g_ub should be filled in here
    // they all have size m
    g_lb.setConstant(-1);
    g_ub.setConstant(1);
}

}; // namespace IDTO

