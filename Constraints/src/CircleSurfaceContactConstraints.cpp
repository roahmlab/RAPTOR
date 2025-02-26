#include "CircleSurfaceContactConstraints.h"

namespace RAPTOR {

CircleSurfaceContactConstraints::CircleSurfaceContactConstraints(std::shared_ptr<CustomizedInverseDynamics>& idPtr_input,
                                                                 const circleContactSurfaceParams& csp_input) :
    idPtr_(idPtr_input),
    csp(csp_input) {
    // (1) positive contact force
    // (2) translation friction cone
    // (3) ZMP (object not flipping over)
    initialize_memory(idPtr_->trajPtr_->N * 3, 
                      idPtr_->trajPtr_->varLength,
                      false);
}

void CircleSurfaceContactConstraints::compute(const VecX& z, 
                                 bool compute_derivatives,
                                 bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    if (compute_hessian) {
        throw std::invalid_argument("CircleSurfaceContactConstraints does not support hessian computation");
    }
    
    idPtr_->compute(z, compute_derivatives, compute_hessian);

    for (int i = 0; i < idPtr_->trajPtr_->N; i++) {
        // access the contact wrench, 
        // which should be the contact wrench between the end effector and the object.
        const Vec6& lambda = idPtr_->lambda(i);

        const Vec3& rotation_torque = lambda.head(3);
        const Vec3& translation_force = lambda.tail(3);

        const double contact_force      = translation_force(2) + csp.maxSuctionForce;
        const double friction_force     = std::sqrt(pow(translation_force(0), 2) + pow(translation_force(1), 2));
        const double ZMP_moment         = std::sqrt(pow(rotation_torque(0), 2) + pow(rotation_torque(1), 2));

        const double max_friction_force = csp.mu * contact_force;
        const double max_ZMP_moment     = csp.R * contact_force;

        // (1) positive contact force
        g(i * 3 + 0) = contact_force;

        // (2) translation friction cone
        g(i * 3 + 1) = friction_force - max_friction_force;

        // (3) ZMP (object not flipping over)
        g(i * 3 + 2) = ZMP_moment - max_ZMP_moment;

        if (compute_derivatives) {
            // assume the contact wrench is always located at the end
            const MatX& protation_torque_pz = idPtr_->plambda_pz(i).topRows(3);
            const MatX& ptranslation_force_pz = idPtr_->plambda_pz(i).bottomRows(3);

            // (1) positive contact force
            pg_pz.row(i * 3 + 0) = ptranslation_force_pz.row(2);

            // (2) translation friction cone
            if (friction_force <= 1e-8) { // avoid singularity when friction_force is close to 0
                pg_pz.row(i * 3 + 1) = -csp.mu * ptranslation_force_pz.row(2);
            }
            else {
                pg_pz.row(i * 3 + 1) = (translation_force(0) * ptranslation_force_pz.row(0) + 
                                        translation_force(1) * ptranslation_force_pz.row(1)) / friction_force - 
                                       csp.mu * ptranslation_force_pz.row(2);
            }

            // (3) ZMP (object not flipping over)
            if (ZMP_moment <= 1e-8) { // avoid singularity when ZMP_moment is close to 0
                pg_pz.row(i * 3 + 2) = -csp.R * ptranslation_force_pz.row(2);
            }
            else {
                pg_pz.row(i * 3 + 2) = (rotation_torque(0) * protation_torque_pz.row(0) + 
                                        rotation_torque(1) * protation_torque_pz.row(1)) / ZMP_moment - 
                                       csp.R * ptranslation_force_pz.row(2);
            }
        }
    }
}

void CircleSurfaceContactConstraints::compute_bounds() {
    for (int i = 0; i < idPtr_->trajPtr_->N; i++) {
        // (1) positive contact force
        g_lb(i * 3 + 0) = csp.contactForceBuffer;
        g_ub(i * 3 + 0) = 1e19;

        // (2) translation friction cone
        g_lb(i * 3 + 1) = -1e19;
        g_ub(i * 3 + 1) = -csp.frictionForceBuffer;

        // (3) ZMP (object not flipping over)
        g_lb(i * 3 + 2) = -1e19; 
        g_ub(i * 3 + 2) = -csp.ZMPBuffer;
    }
}

void CircleSurfaceContactConstraints::print_violation_info() {
    for (int i = 0; i < idPtr_->trajPtr_->N; i++) {
        // (1) positive contact force
        if (g(i * 3 + 0) <= 0) {
            std::cout << "CircleSurfaceContactConstraints: positive contact force violation " 
                      << "at time step " << i << ": " << g(i * 3 + 0) << std::endl;
        }

        // (2) translation friction cone
        if (g(i * 3 + 1) >= 0) {
            std::cout << "CircleSurfaceContactConstraints: translation friction cone violation "
                      << "at time step " << i << ": " << g(i * 3 + 1) << std::endl;
        }

        // (3) ZMP (object not flipping over)
        if (g(i * 3 + 2) >= 0) {
            std::cout << "CircleSurfaceContactConstraints: ZMP (object not flipping over) violation "
                      << "at time step " << i << ": " << g(i * 3 + 2) << std::endl;
        }
    }
}

}; // namespace RAPTOR

