#include "SurfaceContactConstraints.h"

namespace IDTO {

SurfaceContactConstraints::SurfaceContactConstraints(std::shared_ptr<ConstrainedInverseDynamics>& idPtr_input,
                                                     const double friction_mu_input,
                                                     const double friction_gamma_input, 
                                                     const double Lx_input,
                                                     const double Ly_input) : 
    idPtr_(idPtr_input), 
    friction_mu(friction_mu_input),
    friction_gamma(friction_gamma_input),
    Lx(Lx_input),
    Ly(Ly_input) {
    // (1)    positive contact force
    // (2)    translation friction cone
    // (3)    rotation friction cone
    // (4, 5) ZMP on one axis
    // (6, 7) ZMP on the other axis
    m = idPtr_->trajPtr_->N * 7;

    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, idPtr_->trajPtr_->varLength);
}

void SurfaceContactConstraints::compute(const VecX& z, bool compute_derivatives) {
    if (compute_derivatives) {
        pg_pz.setZero();
    }

    // We assume that surface contact constraints always come after torque limits constraints for now
    // The following line has been called in TorqueLimits::compute() already
    // So we directly pull out the lambda values from idPtr_
    // idPtr_->compute(z, compute_derivatives);

    for (int i = 0; i < idPtr_->trajPtr_->N; i++) {
        // assume the contact wrench is always located at the end
        VecX lambda = idPtr_->lambda(i).tail(6);

        double contact_force         = lambda(2);
        double friction_force_sq     = pow(lambda(0), 2) + pow(lambda(1), 2);
        double max_friction_force_sq = pow(friction_mu * contact_force, 2);
        double max_moment_z_sq       = pow(friction_gamma * contact_force, 2);
        double mx_lower_limit        = -Lx * contact_force;
        double mx_upper_limit        = Lx * contact_force;
        double my_lower_limit        = -Ly * contact_force;
        double my_upper_limit        = Ly * contact_force;

        // (1) positive contact force
        g(i * 7 + 0) = contact_force;

        // (2) translation friction cone
        g(i * 7 + 1) = friction_force_sq - max_friction_force_sq;

        // (3) rotation friction cone
        g(i * 7 + 2) = pow(lambda(5), 2) - max_moment_z_sq; 

        // (4, 5) ZMP on one axis
        g(i * 7 + 3) = lambda(3) - mx_upper_limit;
        g(i * 7 + 4) = mx_lower_limit - lambda(3);

        // (6, 7) ZMP on the other axis
        g(i * 7 + 5) = lambda(4) - my_upper_limit;
        g(i * 7 + 6) = my_lower_limit - lambda(4);

        if (compute_derivatives) {
            // assume the contact wrench is always located at the end
            MatX plambda_pz = idPtr_->plambda_pz(i).bottomRows(6);

            // (1) positive contact force
            pg_pz.row(i * 7 + 0) = plambda_pz.row(2);

            // (2) translation friction cone
            pg_pz.row(i * 7 + 1) = 2 * lambda(0) * plambda_pz.row(0) + 
                                   2 * lambda(1) * plambda_pz.row(1) - 
                                   2 * pow(friction_mu, 2) * contact_force * plambda_pz.row(2);

            // (3) rotation friction cone
            pg_pz.row(i * 7 + 2) = 2 * lambda(5) * plambda_pz.row(5) - 
                                   2 * pow(friction_gamma, 2) * contact_force * plambda_pz.row(2);; 

            // (4, 5) ZMP on one axis
            pg_pz.row(i * 7 + 3) = plambda_pz.row(3) - Lx * plambda_pz.row(2);
            pg_pz.row(i * 7 + 4) = -Lx * plambda_pz.row(2) - plambda_pz.row(3);

            // (6, 7) ZMP on the other axis
            pg_pz.row(i * 7 + 5) = plambda_pz.row(4) - Ly * plambda_pz.row(2);
            pg_pz.row(i * 7 + 6) = -Ly * plambda_pz.row(2) - plambda_pz.row(4);
        }
    }
}

void SurfaceContactConstraints::compute_bounds() {
    for (int i = 0; i < idPtr_->trajPtr_->N; i++) {
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

}; // namespace IDTO