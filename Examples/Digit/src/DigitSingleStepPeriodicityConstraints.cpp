#include "DigitSingleStepPeriodicityConstraints.h"

namespace IDTO {
namespace Digit {

DigitSingleStepPeriodicityConstraints::DigitSingleStepPeriodicityConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                                                             std::shared_ptr<DigitConstrainedInverseDynamics> dcidPtr_input) : 
    trajPtr_(trajPtr_input),
    dcidPtr_(dcidPtr_input) {
    m = NUM_JOINTS +              // H * (v+ - v-) = J * lambda 
        NUM_DEPENDENT_JOINTS +    // J*v+ = 0
        NUM_INDEPENDENT_JOINTS +  // position reset
        NUM_INDEPENDENT_JOINTS +  // velocity reset
        7;                        // lambda contact constraints

    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, trajPtr_->varLength);
}

void DigitSingleStepPeriodicityConstraints::compute(const VecX& z, bool compute_derivatives) {
    // We assume that surface contact constraints always come after torque limits constraints for now
    // The following line has been called in TorqueLimits::compute() already
    // So we directly pull out the lambda values from idPtr_
    // dcidPtr_->compute(z, compute_derivatives);

    const VecX& q_minus = dcidPtr_->trajPtr_->q.rightCols(1);
    const VecX& v_minus = dcidPtr_->trajPtr_->q_d.rightCols(1);

    const VecX& q_plus = dcidPtr_->trajPtr_->q.leftCols(1);
    const VecX& v_plus = z.block(trajPtr_->varLength, 0, NUM_JOINTS, 1);

    const VecX& lambda = z.block(trajPtr_->varLength + NUM_JOINTS, 0, NUM_DEPENDENT_JOINTS, 1);

    // swap stance leg for reset map
    dcidPtr_->dcPtr_->reinitialize();

    // re-evaluate constraint jacobian
    dcidPtr_->dcPtr_->get_c(q_minus);
    dcidPtr_->dcPtr_->get_J(q_minus);

    // (1) H * (v+ - v-) = J * lambda
    crba(*(dcidPtr_->modelPtr_), *(dcidPtr_->dataPtr_), q_minus);

    MatX H = dcidPtr_->dataPtr_->M;
    for (size_t j = 0; j < dcidPtr_->modelPtr_->nv; j++) {
        for (size_t k = j + 1; k < dcidPtr_->modelPtr_->nv; k++) {
            H(k, j) = H(j, k);
        }
    }

    VecX g1 = H * (v_plus - v_minus) - dcidPtr_->dcPtr_->J * lambda;

    // (2) J * v+ = 0
    VecX g2 = dcidPtr_->dcPtr_->J * v_plus;

    // (3) position reset
    VecX g3;

    // (4) velocity reset
    VecX g4;

    // (5) lambda contact constraints
    VecX g5(7);
    VecX contactLambda = lambda.tail(6);

    double contact_force         = contactLambda(2);
    double friction_force_sq     = pow(contactLambda(0), 2) + pow(contactLambda(1), 2);
    double max_friction_force_sq = pow(friction_mu * contact_force, 2);
    double max_moment_z_sq       = pow(friction_gamma * contact_force, 2);
    double mx_lower_limit        = -Lx * contact_force;
    double mx_upper_limit        = Lx * contact_force;
    double my_lower_limit        = -Ly * contact_force;
    double my_upper_limit        = Ly * contact_force;

    // (1) positive contact force
    g5(0) = contact_force;

    // (2) translation friction cone
    g5(1) = friction_force_sq - max_friction_force_sq;

    // (3) rotation friction cone
    g5(2) = pow(contactLambda(5), 2) - max_moment_z_sq; 

    // (4, 5) ZMP on one axis
    g5(3) = contactLambda(3) - mx_upper_limit;
    g5(4) = mx_lower_limit - contactLambda(3);

    // (6, 7) ZMP on the other axis
    g5(5) = contactLambda(4) - my_upper_limit;
    g5(6) = my_lower_limit - contactLambda(4);

    // swap back for next round evaluation
    dcidPtr_->dcPtr_->reinitialize();

    // combine all constraints together
    g << g1, g2, g3, g4, g5;
}

void DigitSingleStepPeriodicityConstraints::compute_bounds() {

}

}; // namespace Digit
}; // namespace IDTO