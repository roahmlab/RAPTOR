#include "SurfaceContactConstraints.h"

namespace RAPTOR {

SurfaceContactConstraints::SurfaceContactConstraints(std::shared_ptr<ConstrainedInverseDynamics>& cidPtr_input,
                                                     const frictionParams& fp_input) : 
    cidPtr_(cidPtr_input), 
    fp(fp_input) {
    // (1)    positive contact force
    // (2)    translation friction cone
    // (3)    rotation friction cone
    // (4, 5) ZMP on one axis
    // (6, 7) ZMP on the other axis
    m = cidPtr_->trajPtr_->N * 7;

    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, cidPtr_->trajPtr_->varLength);
}

void SurfaceContactConstraints::compute(const VecX& z, 
                                        bool compute_derivatives,
                                        bool compute_hessian) {
    if (is_computed(z, compute_derivatives, compute_hessian)) {
        return;
    }

    if (compute_hessian) {
        throw std::invalid_argument("SurfaceContactConstraints does not support hessian computation");
    }
    
    cidPtr_->compute(z, compute_derivatives);

    for (int i = 0; i < cidPtr_->trajPtr_->N; i++) {
        // assume the contact wrench is always located at the end
        const VecX& lambda = cidPtr_->lambda(i).tail(6);

        double contact_force      = lambda(2);
        double friction_force     = sqrt(pow(lambda(0), 2) + pow(lambda(1), 2));
        double max_friction_force = fp.mu * contact_force;
        double max_moment_z       = fp.gamma * contact_force;
        double mx_lower_limit     = -fp.Lx * contact_force;
        double mx_upper_limit     = fp.Lx * contact_force;
        double my_lower_limit     = -fp.Ly * contact_force;
        double my_upper_limit     = fp.Ly * contact_force;

        // (1) positive contact force
        g(i * 7 + 0) = contact_force;

        // (2) translation friction cone
        g(i * 7 + 1) = friction_force - max_friction_force;

        // (3) rotation friction cone
        g(i * 7 + 2) = abs(lambda(5)) - max_moment_z; 

        // (4, 5) ZMP on one axis
        g(i * 7 + 3) = lambda(3) - mx_upper_limit;
        g(i * 7 + 4) = mx_lower_limit - lambda(3);

        // (6, 7) ZMP on the other axis
        g(i * 7 + 5) = lambda(4) - my_upper_limit;
        g(i * 7 + 6) = my_lower_limit - lambda(4);

        if (compute_derivatives) {
            // assume the contact wrench is always located at the end
            const MatX& plambda_pz = cidPtr_->plambda_pz(i).bottomRows(6);

            // (1) positive contact force
            pg_pz.row(i * 7 + 0) = plambda_pz.row(2);

            // (2) translation friction cone
            if (friction_force <= 1e-8) {
                pg_pz.row(i * 7 + 1) = fp.mu * plambda_pz.row(2);
            }
            else {
                pg_pz.row(i * 7 + 1) = (lambda(0) * plambda_pz.row(0) + 
                                        lambda(1) * plambda_pz.row(1)) / friction_force - 
                                       fp.mu * plambda_pz.row(2);
            }

            // (3) rotation friction cone
            pg_pz.row(i * 7 + 2) = Utils::sign(lambda(5)) * plambda_pz.row(5) - 
                                   fp.gamma * plambda_pz.row(2);

            // (4, 5) ZMP on one axis
            pg_pz.row(i * 7 + 3) = plambda_pz.row(3) - fp.Lx * plambda_pz.row(2);
            pg_pz.row(i * 7 + 4) = -fp.Lx * plambda_pz.row(2) - plambda_pz.row(3);

            // (6, 7) ZMP on the other axis
            pg_pz.row(i * 7 + 5) = plambda_pz.row(4) - fp.Ly * plambda_pz.row(2);
            pg_pz.row(i * 7 + 6) = -fp.Ly * plambda_pz.row(2) - plambda_pz.row(4);
        }
    }
}

void SurfaceContactConstraints::compute_bounds() {
    for (int i = 0; i < cidPtr_->trajPtr_->N; i++) {
        // (1) positive contact force
        g_lb(i * 7 + 0) = 0;
        g_ub(i * 7 + 0) = 1e19;

        // (2) translation friction cone
        g_lb(i * 7 + 1) = -1e19;
        g_ub(i * 7 + 1) = 0;

        // (3) rotation friction cone
        g_lb(i * 7 + 2) = -1e19; 
        g_ub(i * 7 + 2) = 0; 

        // (4, 5) ZMP on one axis
        g_lb(i * 7 + 3) = -1e19; 
        g_ub(i * 7 + 3) = 0;

        g_lb(i * 7 + 4) = -1e19; 
        g_ub(i * 7 + 4) = 0;

        // (6, 7) ZMP on the other axis
        g_lb(i * 7 + 5) = -1e19; 
        g_ub(i * 7 + 5) = 0;

        g_lb(i * 7 + 6) = -1e19; 
        g_ub(i * 7 + 6) = 0;
    }
}

}; // namespace RAPTOR