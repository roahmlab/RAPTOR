#include "WaitrContactConstraints.h"

namespace RAPTOR {

WaitrContactConstraints::WaitrContactConstraints(std::shared_ptr<CustomizedInverseDynamics>& idPtr_input,
                                                 const contactSurfaceParams& csp_input) :
    idPtr_(idPtr_input),
    csp(csp_input) {
    // (1)    positive contact force
    // (2)    translation friction cone
    // (3, 4) ZMP on one axis
    // (5, 6) ZMP on the other axis
    m = idPtr_->trajPtr_->N * 6; // TODO: this is the number of constraints for all time instances

    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, idPtr_->trajPtr_->varLength);
}

void WaitrContactConstraints::compute(const VecX& z, 
                                      bool compute_derivatives,
                                      bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    if (compute_hessian) {
        throw std::invalid_argument("WaitrContactConstraints does not support hessian computation");
    }
    
    idPtr_->compute(z, compute_derivatives);

    for (int i = 0; i < idPtr_->trajPtr_->N; i++) {
        // access the last contact wrench, 
        // which should be the contact wrench between the end effector and the object
        const Vec6& lambda = idPtr_->lambda(i);

        const Vec3& rotation_torque = lambda.head(3);
        const Vec3& translation_force = lambda.tail(3);

        double contact_force      = translation_force(2);
        double friction_force     = sqrt(pow(translation_force(0), 2) + pow(translation_force(1), 2));
        double max_friction_force = csp.mu * contact_force;
        double mx_lower_limit     = -csp.Lx * contact_force;
        double mx_upper_limit     = csp.Lx * contact_force;
        double my_lower_limit     = -csp.Ly * contact_force;
        double my_upper_limit     = csp.Ly * contact_force;

        // (1) positive contact force
        g(i * 6 + 0) = contact_force;

        // (2) translation friction cone
        g(i * 6 + 1) = friction_force - max_friction_force;

        // (3, 4) ZMP on one axis
        g(i * 6 + 2) = rotation_torque(0) - mx_upper_limit;
        g(i * 6 + 3) = mx_lower_limit - rotation_torque(0);

        // (5, 6) ZMP on the other axis
        g(i * 6 + 4) = rotation_torque(1) - my_upper_limit;
        g(i * 6 + 5) = my_lower_limit - rotation_torque(1);

        if (compute_derivatives) {
            // assume the contact wrench is always located at the end
            const MatX& protation_torque_pz = idPtr_->plambda_pz(i).topRows(3);
            const MatX& ptranslation_force_pz = idPtr_->plambda_pz(i).bottomRows(3);

            // (1) positive contact force
            pg_pz.row(i * 6 + 0) = ptranslation_force_pz.row(2);

            // (2) translation friction cone
            if (friction_force <= 1e-3) { // avoid singularity when friction_force is close to 0
                pg_pz.row(i * 6 + 1) = csp.mu * ptranslation_force_pz.row(2);
            }
            else {
                pg_pz.row(i * 6 + 1) = (translation_force(0) * ptranslation_force_pz.row(0) + 
                                        translation_force(1) * ptranslation_force_pz.row(1)) / friction_force - 
                                       csp.mu * ptranslation_force_pz.row(2);
            }

            // // (3, 4) ZMP on one axis
            pg_pz.row(i * 6 + 2) = protation_torque_pz.row(0) - csp.Lx * ptranslation_force_pz.row(2);
            pg_pz.row(i * 6 + 3) = -csp.Lx * ptranslation_force_pz.row(2) - protation_torque_pz.row(0);

            // (5, 6) ZMP on the other axis
            pg_pz.row(i * 6 + 4) = protation_torque_pz.row(1) - csp.Ly * ptranslation_force_pz.row(2);
            pg_pz.row(i * 6 + 5) = -csp.Ly * ptranslation_force_pz.row(2) - protation_torque_pz.row(1);
        }
    }
}

void WaitrContactConstraints::compute_bounds() {
    for (int i = 0; i < idPtr_->trajPtr_->N; i++) {
        // (1) positive contact force
        g_lb(i * 6 + 0) = 0;
        g_ub(i * 6 + 0) = 1e19;

        // (2) translation friction cone
        g_lb(i * 6 + 1) = -1e19;
        g_ub(i * 6 + 1) = 0;

        // (3, 4) ZMP on one axis
        g_lb(i * 6 + 2) = -1e19; 
        g_ub(i * 6 + 2) = 0;

        g_lb(i * 6 + 3) = -1e19; 
        g_ub(i * 6 + 3) = 0;

        // (5, 6) ZMP on the other axis
        g_lb(i * 6 + 4) = -1e19; 
        g_ub(i * 6 + 4) = 0;

        g_lb(i * 6 + 5) = -1e19; 
        g_ub(i * 6 + 5) = 0;
    }
}

void WaitrContactConstraints::print_violation_info() {
    for (int i = 0; i < idPtr_->trajPtr_->N; i++) {
        if (g(i * 6 + 0) <= 0) {
            std::cout << "        WaitrContactConstraints.cpp: Negative contact force at time instance "
                      << i
                      << std::endl;
        }
        if (g(i * 6 + 1) >= 0) {
            std::cout << "        WaitrContactConstraints.cpp: Translational friction out of friction cone at time instance "
                      << i
                      << std::endl;
        }
        if (g(i * 6 + 2) >= 0) {
            std::cout << "        WaitrContactConstraints.cpp: Object peel off over x dimension at time instance "
                      << i
                      << std::endl;
        }
        if (g(i * 6 + 3) >= 0) {
            std::cout << "        WaitrContactConstraints.cpp: Object peel off over x dimension at time instance "
                      << i
                      << std::endl;
        }
        if (g(i * 6 + 4) >= 0) {
            std::cout << "        WaitrContactConstraints.cpp: Object peel off over y dimension at time instance "
                      << i
                      << std::endl;
        }
        if (g(i * 6 + 5) >= 0) {
            std::cout << "        WaitrContactConstraints.cpp: Object peel off over y dimension at time instance "
                      << i
                      << std::endl;
        }
    }
}

}; // namespace RAPTOR

