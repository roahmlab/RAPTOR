#include "G1MultipleStepPeriodicityConstraints.h"

namespace RAPTOR {
namespace G1 {

G1MultipleStepPeriodicityConstraints::G1MultipleStepPeriodicityConstraints(std::shared_ptr<Trajectories>& currTrajPtr_input,
                                                                           std::shared_ptr<Trajectories>& nextTrajPtr_input,
                                                                           std::shared_ptr<G1ConstrainedInverseDynamics> currDcidPtr_input,
                                                                           std::shared_ptr<G1ConstrainedInverseDynamics> nextDcidPtr_input,
                                                                           const rectangleContactSurfaceParams& fp_input) : 
    currTrajPtr_(currTrajPtr_input),
    nextTrajPtr_(nextTrajPtr_input),
    currDcidPtr_(currDcidPtr_input),
    nextDcidPtr_(nextDcidPtr_input),
    fp(fp_input) {
    if (currTrajPtr_->varLength != nextTrajPtr_->varLength) {
        throw std::invalid_argument("G1MultipleStepPeriodicityConstraints: currTrajPtr_ and nextTrajPtr_ should have the same varLength");
    }

    m = NUM_JOINTS +              // H * (v+ - v-) = J * lambda 
        NUM_DEPENDENT_JOINTS +    // J*v+ = 0
        NUM_INDEPENDENT_JOINTS +  // position reset
        NUM_INDEPENDENT_JOINTS +  // velocity reset
        7;                        // lambda contact constraints
    
    initialize_memory(m, currTrajPtr_->varLength, false);

    // initialize intermediate variables
    const Model& model = *(currDcidPtr_->modelPtr_);
    prnea_pq = MatX::Zero(model.nv, model.nv);
    prnea_pv = MatX::Zero(model.nv, model.nv);
    prnea_pa = MatX::Zero(model.nv, model.nv);
    zeroVec = VecX::Zero(model.nv);

    g1 = VecX::Zero(NUM_JOINTS);
    g2 = VecX::Zero(NUM_DEPENDENT_JOINTS);
    g3 = VecX::Zero(NUM_INDEPENDENT_JOINTS);
    g4 = VecX::Zero(NUM_INDEPENDENT_JOINTS);
    g5 = VecX::Zero(7);

    pg1_pz = MatX::Zero(NUM_JOINTS, currTrajPtr_->varLength);
    pg2_pz = MatX::Zero(NUM_DEPENDENT_JOINTS, currTrajPtr_->varLength);
    pg3_pz = MatX::Zero(NUM_INDEPENDENT_JOINTS, currTrajPtr_->varLength);
    pg3_pz2 = MatX::Zero(NUM_INDEPENDENT_JOINTS, nextTrajPtr_->varLength);
    pg4_pz = MatX::Zero(NUM_INDEPENDENT_JOINTS, currTrajPtr_->varLength);
    pg4_pz2 = MatX::Zero(NUM_INDEPENDENT_JOINTS, nextTrajPtr_->varLength);
    pg5_pz = MatX::Zero(7, currTrajPtr_->varLength);

    pv_plus_pz = MatX::Zero(NUM_JOINTS, currTrajPtr_->varLength);
    for (int i = 0; i < NUM_JOINTS; i++) {
        pv_plus_pz(i, currTrajPtr_->varLength - NUM_JOINTS - NUM_DEPENDENT_JOINTS + i) = 1;
    }

    plambda_pz = MatX::Zero(NUM_DEPENDENT_JOINTS, currTrajPtr_->varLength);
    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        plambda_pz(i, currTrajPtr_->varLength - NUM_DEPENDENT_JOINTS + i) = 1;
    }
}

void G1MultipleStepPeriodicityConstraints::compute(const VecX& z, 
                                                   bool compute_derivatives,
                                                   bool compute_hessian) {
    if (compute_hessian) {
        throw std::invalid_argument("G1MultipleStepPeriodicityConstraints does not support hessian computation");
    }

    Model& model = *(currDcidPtr_->modelPtr_);
    Data& data = *(currDcidPtr_->dataPtr_);
    std::shared_ptr<G1DynamicsConstraints> nextDdcPtr_ = nextDcidPtr_->ddcPtr_;

    const int lastIdx = currTrajPtr_->N - 1;

    const VecX& q_minus = currDcidPtr_->q(lastIdx);
    const MatX& pq_minus_pz = currDcidPtr_->pq_pz(lastIdx);
    
    const VecX& v_minus = currDcidPtr_->v(lastIdx);
    const MatX& pv_minus_pz = currDcidPtr_->pv_pz(lastIdx);

    const VecX& v_plus = z.segment(currTrajPtr_->varLength - NUM_JOINTS - NUM_DEPENDENT_JOINTS, NUM_JOINTS);
    // MatX pv_plus_pz is constant and has been defined in the constructor

    const VecX& q_next = nextDcidPtr_->q(0);
    const MatX& pq_next_pz2 = nextDcidPtr_->pq_pz(0);

    const VecX& v_next = nextDcidPtr_->v(0);
    const MatX& pv_next_pz2 = nextDcidPtr_->pv_pz(0);

    const VecX& lambda = z.tail(NUM_DEPENDENT_JOINTS);

    // evaluate constraint jacobian using next step
    nextDdcPtr_->get_c(q_minus);
    nextDdcPtr_->get_J(q_minus);

    // (1) H * (v+ - v-) = J * lambda
    crba(model, data, q_minus);

    MatX H = data.M;
    for (size_t j = 0; j < model.nv; j++) {
        for (size_t k = j + 1; k < model.nv; k++) {
            H(k, j) = H(j, k);
        }
    }

    g1 = H * (v_plus - v_minus) - nextDdcPtr_->J.transpose() * lambda;

    if (compute_derivatives) {
        // compute derivatives with gravity turned off, we just need prnea_pq here
        const double original_gravity = model.gravity.linear()(2);
        model.gravity.linear()(2) = 0;
        pinocchio::computeRNEADerivatives(model, data, 
                                          q_minus, zeroVec, v_plus - v_minus,
                                          prnea_pq, prnea_pv, prnea_pa);

        // restore gravity
        model.gravity.linear()(2) = original_gravity;

        nextDdcPtr_->get_JTx_partial_dq(q_minus, lambda);
        const MatX& JTx_partial_dq = nextDdcPtr_->JTx_partial_dq;

        pg1_pz = (prnea_pq - JTx_partial_dq) * pq_minus_pz + 
                 H * (pv_plus_pz - pv_minus_pz) - 
                 nextDdcPtr_->J.transpose() * plambda_pz;
    }

    // (2) J * v+ = 0
    g2 = nextDdcPtr_->J * v_plus;

    if (compute_derivatives) {
        nextDdcPtr_->get_Jx_partial_dq(q_minus, v_plus);
        const MatX& Jx_partial_dq = nextDdcPtr_->Jx_partial_dq;

        pg2_pz = Jx_partial_dq * pq_minus_pz + 
                 nextDdcPtr_->J * pv_plus_pz;
    }

    // (3) position reset
    g3 = nextDdcPtr_->get_independent_vector(q_next) - 
         nextDdcPtr_->get_independent_vector(q_minus);

    if (compute_derivatives) {
        nextDdcPtr_->get_independent_rows(pg3_pz2, pq_next_pz2);
        nextDdcPtr_->get_independent_rows(pg3_pz, pq_minus_pz);
        pg3_pz = -pg3_pz;
    }

    // (4) velocity reset
    g4 = nextDdcPtr_->get_independent_vector(v_next) - 
         nextDdcPtr_->get_independent_vector(v_plus);

    if (compute_derivatives) {
        nextDdcPtr_->get_independent_rows(pg4_pz2, pv_next_pz2);
        nextDdcPtr_->get_independent_rows(pg4_pz, pv_plus_pz);
        pg4_pz = -pg4_pz;
    }

    // (5) contact constraints
    VecX contactLambda = lambda.tail(6);

    double contact_force         = contactLambda(2);
    double friction_force_sq     = pow(contactLambda(0), 2) + pow(contactLambda(1), 2);
    double max_friction_force_sq = pow(fp.mu * contact_force, 2);
    double max_moment_z_sq       = pow(fp.gamma * contact_force, 2);
    double mx_lower_limit        = -fp.Lx * contact_force;
    double mx_upper_limit        = fp.Lx * contact_force;
    double my_lower_limit        = -fp.Ly * contact_force;
    double my_upper_limit        = fp.Ly * contact_force;

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

    if (compute_derivatives) {
        // assume the contact wrench is always located at the end
        MatX pcontactLambda_pz = plambda_pz.bottomRows(6);

        // (1) positive contact force
        pg5_pz.row(0) = pcontactLambda_pz.row(2);

        // (2) translation friction cone
        pg5_pz.row(1) = 2 * contactLambda(0) * pcontactLambda_pz.row(0) + 
                        2 * contactLambda(1) * pcontactLambda_pz.row(1) - 
                        2 * pow(fp.mu, 2) * contact_force * pcontactLambda_pz.row(2);

        // (3) rotation friction cone
        pg5_pz.row(2) = 2 * contactLambda(5) * pcontactLambda_pz.row(5) - 
                        2 * pow(fp.gamma, 2) * contact_force * pcontactLambda_pz.row(2);; 

        // (4, 5) ZMP on one axis
        pg5_pz.row(3) = pcontactLambda_pz.row(3) - fp.Lx * pcontactLambda_pz.row(2);
        pg5_pz.row(4) = -fp.Lx * pcontactLambda_pz.row(2) - pcontactLambda_pz.row(3);

        // (6, 7) ZMP on the other axis
        pg5_pz.row(5) = pcontactLambda_pz.row(4) - fp.Ly * pcontactLambda_pz.row(2);
        pg5_pz.row(6) = -fp.Ly * pcontactLambda_pz.row(2) - pcontactLambda_pz.row(4);
    }

    // combine all constraints together
    g << g1, g2, g3, g4, g5;

    if (compute_derivatives) {
        pg_pz << pg1_pz, pg2_pz, pg3_pz, pg4_pz, pg5_pz;

        // pg3_pz2, pg4_pz2 is for the next walking step
        // they are directly extracted from this class and then composed to the values passed to ipopt
    }
}

void G1MultipleStepPeriodicityConstraints::compute_bounds() {
    // everything before this are equality constraints, so just zeros
    // and g_lb, g_ub are already initialized to zeros in the constructor

    // g_lb(NUM_JOINTS + NUM_DEPENDENT_JOINTS + 2 * NUM_INDEPENDENT_JOINTS) = 0;
    g_ub(NUM_JOINTS + NUM_DEPENDENT_JOINTS + 2 * NUM_INDEPENDENT_JOINTS) = 1e19;

    g_lb.segment(NUM_JOINTS + NUM_DEPENDENT_JOINTS + 2 * NUM_INDEPENDENT_JOINTS + 1, 6).setConstant(-1e19);
    // g_ub.segment(NUM_JOINTS + NUM_DEPENDENT_JOINTS + 2 * NUM_INDEPENDENT_JOINTS + 1, 6).setZero();
}

void G1MultipleStepPeriodicityConstraints::print_violation_info() {
    // (1) H * (v+ - v-) = J * lambda
    for (int i = 0; i < NUM_JOINTS; i++) {
        if (abs(g1(i)) >= 1e-4) {
            std::cout << "        G1MultipleStepPeriodicityConstraints.cpp: H * (v+ - v-) = J * lambda: dim " 
                      << i 
                      << " is violated: "
                      << g(i) 
                      << std::endl;
        }
    }

    // (2) J * v+ = 0
    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        if (abs(g2(i)) >= 1e-4) {
            std::cout << "        G1MultipleStepPeriodicityConstraints.cpp: J * v+ = 0: dim " 
                      << i
                      << " is violated: "
                      << g2(i)
                      << std::endl;
        }
    }

    // (3) position reset
    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        if (abs(g3(i)) >= 1e-4) {
            std::cout << "        G1MultipleStepPeriodicityConstraints.cpp: position reset: dim " 
                      << i
                      << " is violated: "
                      << g3(i)
                      << std::endl;
        }
    }

    // (4) velocity reset
    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        if (abs(g4(i)) >= 1e-4) {
            std::cout << "        G1MultipleStepPeriodicityConstraints.cpp: velocity reset: dim " 
                      << i
                      << " is violated: "
                      << g4(i)
                      << std::endl;
        }
    }

    // (5) contact constraints
    if (g5(0) <= -1e-4) {
        std::cout << "        G1MultipleStepPeriodicityConstraints.cpp: positive contact force is violated: " 
                  << g5(0)
                  << std::endl;
    }
    if (g5(1) >= 1e-4) {
        std::cout << "        G1MultipleStepPeriodicityConstraints.cpp: translation friction cone is violated: " 
                  << g5(1)
                  << std::endl;
    }
    if (g5(2) >= 1e-4) {
        std::cout << "        G1MultipleStepPeriodicityConstraints.cpp: rotation friction cone is violated: " 
                  << g5(2)
                  << std::endl;
    }
    if (g5(3) >= 1e-4) {
        std::cout << "        G1MultipleStepPeriodicityConstraints.cpp: ZMP on one axis is violated: " 
                  << g5(3)
                  << std::endl;
    }
    if (g5(4) >= 1e-4) {
        std::cout << "        G1MultipleStepPeriodicityConstraints.cpp: ZMP on one axis is violated: " 
                  << g5(4)
                  << std::endl;
    }
    if (g5(5) >= 1e-4) {
        std::cout << "        G1MultipleStepPeriodicityConstraints.cpp: ZMP on the other axis is violated: " 
                  << g5(5)
                  << std::endl;
    }
    if (g5(6) >= 1e-4) {
        std::cout << "        G1MultipleStepPeriodicityConstraints.cpp: ZMP on the other axis is violated: " 
                  << g5(6)
                  << std::endl;
    }
}

}; // namespace G1
}; // namespace RAPTOR