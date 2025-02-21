#include "BipedSingleStepPeriodicityConstraints.h"

namespace RAPTOR {
namespace Biped {

BipedSingleStepPeriodicityConstraints::BipedSingleStepPeriodicityConstraints(std::shared_ptr<Trajectories>& trajPtr_input,
                                                                             std::shared_ptr<BipedConstrainedInverseDynamics> dcidPtr_input,
                                                                             const rectangleContactSurfaceParams& fp_input) : 
    trajPtr_(trajPtr_input),
    dcidPtr_(dcidPtr_input),
    fp(fp_input) {
    m = NUM_JOINTS +              // H * (v+ - v-) = J * lambda 
        NUM_DEPENDENT_JOINTS +    // J*v+ = 0
        NUM_INDEPENDENT_JOINTS +  // position reset
        NUM_INDEPENDENT_JOINTS +  // velocity reset
        7;                        // lambda contact constraints
    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, trajPtr_->varLength);

    // initialize intermediate variables
    prnea_pq = MatX::Zero(dcidPtr_->modelPtr_->nv, dcidPtr_->modelPtr_->nv);
    prnea_pv = MatX::Zero(dcidPtr_->modelPtr_->nv, dcidPtr_->modelPtr_->nv);
    prnea_pa = MatX::Zero(dcidPtr_->modelPtr_->nv, dcidPtr_->modelPtr_->nv);
    zeroVec = VecX::Zero(dcidPtr_->modelPtr_->nv);

    g1 = VecX::Zero(NUM_JOINTS);
    g2 = VecX::Zero(NUM_DEPENDENT_JOINTS);
    g3 = VecX::Zero(NUM_INDEPENDENT_JOINTS);
    g4 = VecX::Zero(NUM_INDEPENDENT_JOINTS);
    g5 = VecX::Zero(7);

    pg1_pz = MatX::Zero(NUM_JOINTS, trajPtr_->varLength);
    pg2_pz = MatX::Zero(NUM_DEPENDENT_JOINTS, trajPtr_->varLength);
    pg3_pz = MatX::Zero(NUM_INDEPENDENT_JOINTS, trajPtr_->varLength);
    pg4_pz = MatX::Zero(NUM_INDEPENDENT_JOINTS, trajPtr_->varLength);
    pg5_pz = MatX::Zero(7, trajPtr_->varLength);

    // the following are going to be constant
    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        joint_id1[i] = dcidPtr_->modelPtr_->getJointId(JOINT_MAP[i][0]) - 1;
        joint_id2[i] = dcidPtr_->modelPtr_->getJointId(JOINT_MAP[i][1]) - 1;
    }

    pv_plus_pz = MatX::Zero(NUM_JOINTS, trajPtr_->varLength);
    for (int i = 0; i < NUM_JOINTS; i++) {
        pv_plus_pz(i, trajPtr_->varLength - NUM_JOINTS - NUM_DEPENDENT_JOINTS + i) = 1;
    }

    plambda_pz = MatX::Zero(NUM_DEPENDENT_JOINTS, trajPtr_->varLength);
    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        plambda_pz(i, trajPtr_->varLength - NUM_DEPENDENT_JOINTS + i) = 1;
    }
}

void BipedSingleStepPeriodicityConstraints::compute(const VecX& z, 
                                                 bool compute_derivatives,
                                                 bool compute_hessian) {
    if (compute_hessian) {
        throw std::invalid_argument("BipedSingleStepPeriodicityConstraints does not support hessian computation");
    }

    // We assume that surface contact constraints always come after torque limits constraints for now
    // The following line has been called in TorqueLimits::compute() already
    // So we directly pull out the lambda values from idPtr_
    dcidPtr_->compute(z, compute_derivatives);

    const int lastIdx = trajPtr_->N - 1;

    const VecX& q_minus = dcidPtr_->q(lastIdx);
    const VecX& v_minus = dcidPtr_->v(lastIdx);

    const VecX& v_plus = z.segment(trajPtr_->varLength - NUM_JOINTS - NUM_DEPENDENT_JOINTS, NUM_JOINTS);

    const VecX& q_0 = dcidPtr_->q(0);
    const VecX& v_0 = dcidPtr_->v(0);

    const VecX& lambda = z.tail(NUM_DEPENDENT_JOINTS);

    // swap stance leg for reset map
    dcidPtr_->ddcPtr_->reinitialize();

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

    g1 = H * (v_plus - v_minus) - dcidPtr_->dcPtr_->J.transpose() * lambda;

    if (compute_derivatives) {
        // compute derivatives with gravity turned off, we just need prnea_pq here
        const double original_gravity = dcidPtr_->modelPtr_->gravity.linear()(2);
        dcidPtr_->modelPtr_->gravity.linear()(2) = 0;
        pinocchio::computeRNEADerivatives(*(dcidPtr_->modelPtr_), *(dcidPtr_->dataPtr_), 
                                          q_minus, zeroVec, v_plus - v_minus,
                                          prnea_pq, prnea_pv, prnea_pa);

        // restore gravity
        dcidPtr_->modelPtr_->gravity.linear()(2) = original_gravity;

        dcidPtr_->dcPtr_->get_JTx_partial_dq(q_minus, lambda);
        const MatX& JTx_partial_dq = dcidPtr_->dcPtr_->JTx_partial_dq;

        pg1_pz = (prnea_pq - JTx_partial_dq) * dcidPtr_->pq_pz(lastIdx) + 
                 H * (pv_plus_pz - dcidPtr_->pv_pz(lastIdx)) - 
                 dcidPtr_->dcPtr_->J.transpose() * plambda_pz;
    }

    // (2) J * v+ = 0
    g2 = dcidPtr_->dcPtr_->J * v_plus;

    if (compute_derivatives) {
        dcidPtr_->dcPtr_->get_Jx_partial_dq(q_minus, v_plus);
        const MatX& Jx_partial_dq = dcidPtr_->dcPtr_->Jx_partial_dq;

        pg2_pz = Jx_partial_dq * dcidPtr_->pq_pz(lastIdx) + 
                 dcidPtr_->dcPtr_->J * pv_plus_pz;
    }

    // (3) position reset
    // for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
    //     g3(i) = q_0(joint_id1[i]) - q_minus(joint_id2[i]);
    // }
    g3(0) = q_0(joint_id1[0]) + q_minus(joint_id2[0]); // "left_leg_1" ==> "right_leg_1"
    g3(1) = q_0(joint_id1[1]) - q_minus(joint_id2[1]); // "left_leg_2" ==> "right_leg_2"
    g3(2) = q_0(joint_id1[2]) - q_minus(joint_id2[2]); // "left_leg_3" ==> "right_leg_3"
    g3(3) = q_0(joint_id1[3]) + q_minus(joint_id2[3]); // "right_leg_1" ==> "left_leg_1"
    g3(4) = q_0(joint_id1[4]) - q_minus(joint_id2[4]); // "right_leg_2" ==> "left_leg_2"
    g3(5) = q_0(joint_id1[5]) - q_minus(joint_id2[5]); // "right_leg_3" ==> "left_leg_3"
    g3(6) = q_0(joint_id1[6]) + q_minus(joint_id2[6]); // "shoulders" ==> "shoulders"
    g3(7) = q_0(joint_id1[7]) - q_minus(joint_id2[7]); // "head" ==> "head"
    g3(8) = q_0(joint_id1[8]) + q_minus(joint_id2[8]); // "left_arm_1" ==> "right_arm_1"
    g3(9) = q_0(joint_id1[9]) - q_minus(joint_id2[9]); // "left_arm_2" ==> "right_arm_2"
    g3(10) = q_0(joint_id1[10]) - q_minus(joint_id2[10]); // "left_arm_3" ==> "right_arm_3"
    g3(11) = q_0(joint_id1[11]) + q_minus(joint_id2[11]); // "right_arm_1" ==> "left_arm_1"
    g3(12) = q_0(joint_id1[12]) - q_minus(joint_id2[12]); // "right_arm_2" ==> "left_arm_2"
    g3(13) = q_0(joint_id1[13]) - q_minus(joint_id2[13]); // "right_arm_3" ==> "left_arm_3"

    if (compute_derivatives) {
        pg3_pz.row(0) = dcidPtr_->pq_pz(0).row(joint_id1[0]) + dcidPtr_->pq_pz(lastIdx).row(joint_id2[0]);
        pg3_pz.row(1) = dcidPtr_->pq_pz(0).row(joint_id1[1]) - dcidPtr_->pq_pz(lastIdx).row(joint_id2[1]);
        pg3_pz.row(2) = dcidPtr_->pq_pz(0).row(joint_id1[2]) - dcidPtr_->pq_pz(lastIdx).row(joint_id2[2]);
        pg3_pz.row(3) = dcidPtr_->pq_pz(0).row(joint_id1[3]) + dcidPtr_->pq_pz(lastIdx).row(joint_id2[3]);
        pg3_pz.row(4) = dcidPtr_->pq_pz(0).row(joint_id1[4]) - dcidPtr_->pq_pz(lastIdx).row(joint_id2[4]);
        pg3_pz.row(5) = dcidPtr_->pq_pz(0).row(joint_id1[5]) - dcidPtr_->pq_pz(lastIdx).row(joint_id2[5]);
        pg3_pz.row(6) = dcidPtr_->pq_pz(0).row(joint_id1[6]) + dcidPtr_->pq_pz(lastIdx).row(joint_id2[6]);
        pg3_pz.row(7) = dcidPtr_->pq_pz(0).row(joint_id1[7]) - dcidPtr_->pq_pz(lastIdx).row(joint_id2[7]);
        pg3_pz.row(8) = dcidPtr_->pq_pz(0).row(joint_id1[8]) + dcidPtr_->pq_pz(lastIdx).row(joint_id2[8]);
        pg3_pz.row(9) = dcidPtr_->pq_pz(0).row(joint_id1[9]) - dcidPtr_->pq_pz(lastIdx).row(joint_id2[9]);
        pg3_pz.row(10) = dcidPtr_->pq_pz(0).row(joint_id1[10]) - dcidPtr_->pq_pz(lastIdx).row(joint_id2[10]);
        pg3_pz.row(11) = dcidPtr_->pq_pz(0).row(joint_id1[11]) + dcidPtr_->pq_pz(lastIdx).row(joint_id2[11]);
        pg3_pz.row(12) = dcidPtr_->pq_pz(0).row(joint_id1[12]) - dcidPtr_->pq_pz(lastIdx).row(joint_id2[12]);
        pg3_pz.row(13) = dcidPtr_->pq_pz(0).row(joint_id1[13]) - dcidPtr_->pq_pz(lastIdx).row(joint_id2[13]);
    }

    // (4) velocity reset
    // for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
    //     g4(i) = v_0(joint_id1[i]) - v_plus(joint_id2[i]);
    // }
    g4(0) = v_0(joint_id1[0]) + v_plus(joint_id2[0]); // "left_leg_1" ==> "right_leg_1"
    g4(1) = v_0(joint_id1[1]) - v_plus(joint_id2[1]); // "left_leg_2" ==> "right_leg_2"
    g4(2) = v_0(joint_id1[2]) - v_plus(joint_id2[2]); // "left_leg_3" ==> "right_leg_3"
    g4(3) = v_0(joint_id1[3]) + v_plus(joint_id2[3]); // "right_leg_1" ==> "left_leg_1"
    g4(4) = v_0(joint_id1[4]) - v_plus(joint_id2[4]); // "right_leg_2" ==> "left_leg_2"
    g4(5) = v_0(joint_id1[5]) - v_plus(joint_id2[5]); // "right_leg_3" ==> "left_leg_3"
    g4(6) = v_0(joint_id1[6]) + v_plus(joint_id2[6]); // "shoulders" ==> "shoulders"
    g4(7) = v_0(joint_id1[7]) - v_plus(joint_id2[7]); // "head" ==> "head"
    g4(8) = v_0(joint_id1[8]) + v_plus(joint_id2[8]); // "left_arm_1" ==> "right_arm_1"
    g4(9) = v_0(joint_id1[9]) - v_plus(joint_id2[9]); // "left_arm_2" ==> "right_arm_2"
    g4(10) = v_0(joint_id1[10]) - v_plus(joint_id2[10]); // "left_arm_3" ==> "right_arm_3"
    g4(11) = v_0(joint_id1[11]) + v_plus(joint_id2[11]); // "right_arm_1" ==> "left_arm_1"
    g4(12) = v_0(joint_id1[12]) - v_plus(joint_id2[12]); // "right_arm_2" ==> "left_arm_2"
    g4(13) = v_0(joint_id1[13]) - v_plus(joint_id2[13]); // "right_arm_3" ==> "left_arm_3"

    if (compute_derivatives) {
        pg4_pz.row(0) = dcidPtr_->pv_pz(0).row(joint_id1[0]) + pv_plus_pz.row(joint_id2[0]);
        pg4_pz.row(1) = dcidPtr_->pv_pz(0).row(joint_id1[1]) - pv_plus_pz.row(joint_id2[1]);
        pg4_pz.row(2) = dcidPtr_->pv_pz(0).row(joint_id1[2]) - pv_plus_pz.row(joint_id2[2]);
        pg4_pz.row(3) = dcidPtr_->pv_pz(0).row(joint_id1[3]) + pv_plus_pz.row(joint_id2[3]);
        pg4_pz.row(4) = dcidPtr_->pv_pz(0).row(joint_id1[4]) - pv_plus_pz.row(joint_id2[4]);
        pg4_pz.row(5) = dcidPtr_->pv_pz(0).row(joint_id1[5]) - pv_plus_pz.row(joint_id2[5]);
        pg4_pz.row(6) = dcidPtr_->pv_pz(0).row(joint_id1[6]) + pv_plus_pz.row(joint_id2[6]);
        pg4_pz.row(7) = dcidPtr_->pv_pz(0).row(joint_id1[7]) - pv_plus_pz.row(joint_id2[7]);
        pg4_pz.row(8) = dcidPtr_->pv_pz(0).row(joint_id1[8]) + pv_plus_pz.row(joint_id2[8]);
        pg4_pz.row(9) = dcidPtr_->pv_pz(0).row(joint_id1[9]) - pv_plus_pz.row(joint_id2[9]);
        pg4_pz.row(10) = dcidPtr_->pv_pz(0).row(joint_id1[10]) - pv_plus_pz.row(joint_id2[10]);
        pg4_pz.row(11) = dcidPtr_->pv_pz(0).row(joint_id1[11]) + pv_plus_pz.row(joint_id2[11]);
        pg4_pz.row(12) = dcidPtr_->pv_pz(0).row(joint_id1[12]) - pv_plus_pz.row(joint_id2[12]);
        pg4_pz.row(13) = dcidPtr_->pv_pz(0).row(joint_id1[13]) - pv_plus_pz.row(joint_id2[13]);
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

    // swap back for next round evaluation
    dcidPtr_->ddcPtr_->reinitialize();

    // combine all constraints together
    g << g1, g2, g3, g4, g5;

    if (compute_derivatives) {
        pg_pz << pg1_pz, pg2_pz, pg3_pz, pg4_pz, pg5_pz;
    }
}

void BipedSingleStepPeriodicityConstraints::compute_bounds() {
    // everything before this are equality constraints, so just zeros
    // and g_lb, g_ub are already initialized to zeros in the constructor

    // g_lb(NUM_JOINTS + NUM_DEPENDENT_JOINTS + 2 * NUM_INDEPENDENT_JOINTS) = 0;
    g_ub(NUM_JOINTS + NUM_DEPENDENT_JOINTS + 2 * NUM_INDEPENDENT_JOINTS) = 1e19;

    g_lb.segment(NUM_JOINTS + NUM_DEPENDENT_JOINTS + 2 * NUM_INDEPENDENT_JOINTS + 1, 6).setConstant(-1e19);
    // g_ub.segment(NUM_JOINTS + NUM_DEPENDENT_JOINTS + 2 * NUM_INDEPENDENT_JOINTS + 1, 6).setZero();
}

void BipedSingleStepPeriodicityConstraints::print_violation_info() {
    // (1) H * (v+ - v-) = J * lambda
    for (int i = 0; i < NUM_JOINTS; i++) {
        if (abs(g1(i)) > 1e-5) {
            std::cout << "        BipedSingleStepPeriodicityConstraints.cpp: H * (v+ - v-) = J * lambda: dim " 
                      << i 
                      << " is violated: "
                      << g(i) 
                      << std::endl;
        }
    }

    // (2) J * v+ = 0
    for (int i = 0; i < NUM_DEPENDENT_JOINTS; i++) {
        if (abs(g2(i)) > 1e-5) {
            std::cout << "        BipedSingleStepPeriodicityConstraints.cpp: J * v+ = 0: dim " 
                      << i
                      << " is violated: "
                      << g2(i)
                      << std::endl;
        }
    }

    // (3) position reset
    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        if (abs(g3(i)) > 1e-5) {
            std::cout << "        BipedSingleStepPeriodicityConstraints.cpp: position reset: dim " 
                      << i
                      << " is violated: "
                      << g3(i)
                      << std::endl;
        }
    }

    // (4) velocity reset
    for (int i = 0; i < NUM_INDEPENDENT_JOINTS; i++) {
        if (abs(g4(i)) > 1e-5) {
            std::cout << "        BipedSingleStepPeriodicityConstraints.cpp: velocity reset: dim " 
                      << i
                      << " is violated: "
                      << g4(i)
                      << std::endl;
        }
    }

    // (5) contact constraints
    if (g5(0) <= -1e-4) {
        std::cout << "        BipedSingleStepPeriodicityConstraints.cpp: positive contact force is violated: " 
                  << g5(0)
                  << std::endl;
    }
    if (g5(1) >= 1e-4) {
        std::cout << "        BipedSingleStepPeriodicityConstraints.cpp: translation friction cone is violated: " 
                  << g5(1)
                  << std::endl;
    }
    if (g5(2) >= 1e-4) {
        std::cout << "        BipedSingleStepPeriodicityConstraints.cpp: rotation friction cone is violated: " 
                  << g5(2)
                  << std::endl;
    }
    if (g5(3) >= 1e-4) {
        std::cout << "        BipedSingleStepPeriodicityConstraints.cpp: ZMP on one axis is violated: " 
                  << g5(3)
                  << std::endl;
    }
    if (g5(4) >= 1e-4) {
        std::cout << "        BipedSingleStepPeriodicityConstraints.cpp: ZMP on one axis is violated: " 
                  << g5(4)
                  << std::endl;
    }
    if (g5(5) >= 1e-4) {
        std::cout << "        BipedSingleStepPeriodicityConstraints.cpp: ZMP on the other axis is violated: " 
                  << g5(5)
                  << std::endl;
    }
    if (g5(6) >= 1e-4) {
        std::cout << "        BipedSingleStepPeriodicityConstraints.cpp: ZMP on the other axis is violated: " 
                  << g5(6)
                  << std::endl;
    }
}

}; // namespace Biped
}; // namespace RAPTOR