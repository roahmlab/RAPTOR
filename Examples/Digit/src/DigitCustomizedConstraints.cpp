#include "DigitCustomizedConstraints.h"

namespace RAPTOR {
namespace Digit {

DigitCustomizedConstraints::DigitCustomizedConstraints(const Model& model_input,
                                                       std::shared_ptr<Trajectories>& trajPtr_input,
                                                       std::shared_ptr<DigitDynamicsConstraints>& ddcPtr_input,
                                                       const GaitParameters& gp_input) : 
    trajPtr_(trajPtr_input),
    ddcPtr_(ddcPtr_input),
    gp(gp_input) {
    modelPtr_ = std::make_unique<Model>(model_input);
    fkPtr_ = std::make_unique<ForwardKinematicsSolver>(modelPtr_.get());

    leftfoot_endT.R << 0,             1, 0,
                       -0.5,          0, sin(M_PI / 3),
                       sin(M_PI / 3), 0, 0.5;
    leftfoot_endT.p << 0, -0.05456, -0.0315;

    rightfoot_endT.R << 0,             -1, 0,
                        0.5,           0,  -sin(M_PI / 3),
                        sin(M_PI / 3), 0,  0.5;
    rightfoot_endT.p << 0, 0.05456, -0.0315;

    q = MatX::Zero(modelPtr_->nv, trajPtr_->N);
    pq_pz.resize(1, trajPtr_->N);
    swingfoot_xyzrpy = MatX::Zero(6, trajPtr_->N);
    pswingfoot_xyzrpy_pz.resize(1, trajPtr_->N);

    m = trajPtr_->N * 10 + 4;
    initialize_memory(m, trajPtr_->varLength, false);
}

void DigitCustomizedConstraints::compute(const VecX& z, 
                                         bool compute_derivatives,
                                         bool compute_hessian) {
    if (compute_hessian) {
        throw std::invalid_argument("DigitCustomizedConstraints does not support hessian computation");
    }

    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    const int fk_order = compute_derivatives ? 1 : 0;

    // compute full joint trajectories and swing foot forward kinematics
    for (int i = 0; i < trajPtr_->N; i++) {
        if (ddcPtr_->stanceLeg == 'L') {
            swingfoot_endT = rightfoot_endT;
            swingfoot_id = modelPtr_->getJointId(std::string(RIGHT_FOOT_NAME));
        }
        else {
            swingfoot_endT = leftfoot_endT;
            swingfoot_id = modelPtr_->getJointId(std::string(LEFT_FOOT_NAME));
        }
        
        VecX qi(modelPtr_->nq);
        ddcPtr_->fill_independent_vector(qi, trajPtr_->q(i));
        ddcPtr_->setupJointPosition(qi, compute_derivatives);
        q.col(i) = qi;

        fkPtr_->compute(0, swingfoot_id, qi, nullptr, &swingfoot_endT, fk_order);

        swingfoot_xyzrpy.col(i).head(3) = fkPtr_->getTranslation();
        swingfoot_xyzrpy.col(i).tail(3) = fkPtr_->getRPY();

        if (compute_derivatives) {
            pq_pz(i).resize(modelPtr_->nv, trajPtr_->varLength);
            
            // fill in independent joints derivatives directly
            for (int j = 0; j < ddcPtr_->numIndependentJoints; j++) {
                int indenpendentJointIndex = ddcPtr_->return_independent_joint_index(j);
                pq_pz(i).row(indenpendentJointIndex) = trajPtr_->pq_pz(i).row(j);
            }
            // compute and fill in dependent joints derivatives
            MatX pq_dep_pz = ddcPtr_->pq_dep_pq_indep * trajPtr_->pq_pz(i);
            for (int j = 0; j < ddcPtr_->numDependentJoints; j++) {
                int denpendentJointIndex = ddcPtr_->return_dependent_joint_index(j);
                pq_pz(i).row(denpendentJointIndex) = pq_dep_pz.row(j);
            }

            pswingfoot_xyzrpy_pz(i).resize(6, trajPtr_->varLength);
            pswingfoot_xyzrpy_pz(i).topRows(3) = fkPtr_->getTranslationJacobian() * pq_pz(i);
            pswingfoot_xyzrpy_pz(i).bottomRows(3) = fkPtr_->getRPYJacobian() * pq_pz(i);
        }
    }

    // (1) swingfoot height always larger than 0
    //               height equal to 0 at the beginning and at the end
    //               height higher than the desired value in the middle
    g1 = swingfoot_xyzrpy.row(2);

    // (2) swingfoot always flat and points forward
    g2 = Utils::wrapToPi(swingfoot_xyzrpy.row(3)); // swingfoot roll
    g3 = Utils::wrapToPi(swingfoot_xyzrpy.row(4)); // swingfoot pitch
    g4 = (ddcPtr_->stanceLeg == 'L') ? 
        Utils::wrapToPi(swingfoot_xyzrpy.row(5).array() + M_PI_2) :
        Utils::wrapToPi(swingfoot_xyzrpy.row(5).array() - M_PI_2); // swingfoot yaw

    // (3) swingfoot xy equal to desired value 
    //     at the beginning and at the end of each walking step
    g5 = VecX::Zero(4);
    g5 << swingfoot_xyzrpy.topLeftCorner(2, 1), swingfoot_xyzrpy.topRightCorner(2, 1);

    // (4) torso height always larger than 1 meter
    //           torso stays between the two feet
    //           roll and pitch always close to 0
    //           yaw always close to 0 when walking forward
    g6 = q.row(2); // torso height
    g7 = q.row(0) - 0.5 * swingfoot_xyzrpy.row(0); // torso stays between the two feet
    g8 = q.row(1) - 0.5 * swingfoot_xyzrpy.row(1); // torso stays between the two feet
    g9 = Utils::wrapToPi(q.row(3)); // torso roll
    g10 = Utils::wrapToPi(q.row(4)); // torso pitch
    g11 = (ddcPtr_->stanceLeg == 'L') ?
        Utils::wrapToPi(q.row(5).array() + M_PI_2) :
        Utils::wrapToPi(q.row(5).array() - M_PI_2); // torso yaw

    g << g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11;

    if (compute_derivatives) {
        int iter = 0;
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(i).row(2);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(i).row(3);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(i).row(4);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(i).row(5);
        }
        pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(0).row(0);
        pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(0).row(1);
        pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(trajPtr_->N - 1).row(0);
        pg_pz.row(iter++) = pswingfoot_xyzrpy_pz(trajPtr_->N - 1).row(1);
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pq_pz(i).row(2);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pq_pz(i).row(0) - 0.5 * pswingfoot_xyzrpy_pz(i).row(0);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pq_pz(i).row(1) - 0.5 * pswingfoot_xyzrpy_pz(i).row(1);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pq_pz(i).row(3);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pq_pz(i).row(4);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pq_pz(i).row(5);
        }
    }
}

void DigitCustomizedConstraints::compute_bounds() {
    // (1) swingfoot height always larger than 0
    //               height equal to 0 at the beginning and at the end
    //               height higher than the desired value in the middle
    g1_lb = VecX::Zero(trajPtr_->N);
    g1_ub = VecX::Constant(trajPtr_->N, 1e19);
    g1_lb(trajPtr_->N / 2) = gp.swingfoot_midstep_z_des;
    g1_ub(0) = 0;
    g1_ub(trajPtr_->N - 1) = 0;

    // (2) swingfoot always flat and points forward
    g2_lb = VecX::Zero(trajPtr_->N);
    g2_ub = VecX::Zero(trajPtr_->N);
    g3_lb = VecX::Zero(trajPtr_->N);
    g3_ub = VecX::Zero(trajPtr_->N);
    g4_lb = VecX::Zero(trajPtr_->N);
    g4_ub = VecX::Zero(trajPtr_->N);

    // (3) swingfoot xy equal to desired value at the beginning and at the end
    g5_lb = VecX::Zero(4);
    g5_ub = VecX::Zero(4);
    g5_lb << gp.swingfoot_begin_x_des, gp.swingfoot_begin_y_des, gp.swingfoot_end_x_des, gp.swingfoot_end_y_des;
    g5_ub << gp.swingfoot_begin_x_des, gp.swingfoot_begin_y_des, gp.swingfoot_end_x_des, gp.swingfoot_end_y_des;

    // (4) torso height always larger than 1 meter
    //     torso stays between the two feet
    //     roll and pitch always close to 0
    //     yaw always close to 0 when walking forward
    g6_lb = VecX::Constant(trajPtr_->N, 1.05);
    g6_ub = VecX::Constant(trajPtr_->N, 1e19);
    g7_lb = VecX::Constant(trajPtr_->N, -0.05);
    g7_ub = VecX::Constant(trajPtr_->N, 0.05);
    g8_lb = VecX::Constant(trajPtr_->N, -0.10);
    g8_ub = VecX::Constant(trajPtr_->N, 0.10);
    g9_lb = VecX::Constant(trajPtr_->N, -gp.eps_torso_angle);
    g9_ub = VecX::Constant(trajPtr_->N, gp.eps_torso_angle);
    g10_lb = VecX::Constant(trajPtr_->N, -gp.eps_torso_angle);
    g10_ub = VecX::Constant(trajPtr_->N, gp.eps_torso_angle);
    g11_lb = VecX::Constant(trajPtr_->N, -gp.eps_torso_angle);
    g11_ub = VecX::Constant(trajPtr_->N, gp.eps_torso_angle);

    g_lb << g1_lb, g2_lb, g3_lb, g4_lb, g5_lb, g6_lb, g7_lb, g8_lb, g9_lb, g10_lb, g11_lb;
    g_ub << g1_ub, g2_ub, g3_ub, g4_ub, g5_ub, g6_ub, g7_ub, g8_ub, g9_ub, g10_ub, g11_ub;
}

void DigitCustomizedConstraints::print_violation_info() {
    // (1) swingfoot height always larger than 0
    //               height equal to 0 at the beginning and at the end
    //               height higher than the desired value in the middle
    for (int i = 0; i < trajPtr_->N; i++) {
        if (g1(i) <= g1_lb(i) - 1e-4) {
            std::cout << "        DigitCustomizedConstraints.cpp: swing foot height at time instance " 
                      << i 
                      << " is below lower limit: " 
                      << g1(i) 
                      << " < " 
                      << g1_lb(i) 
                      << std::endl;
        }
        if (g1(i) >= g1_ub(i) + 1e-4) {
            std::cout << "        DigitCustomizedConstraints.cpp: swing foot height at time instance " 
                      << i 
                      << " is above upper limit: " 
                      << g1(i) 
                      << " > " 
                      << g1_ub(i) 
                      << std::endl;
        }
    }

    // (2) swingfoot always flat and points forward
    for (int i = 0; i < trajPtr_->N; i++) {
        if (abs(g2(i)) >= 1e-4) {
            std::cout << "        DigitCustomizedConstraints.cpp: swing foot roll at time instance " 
                      << i 
                      << " is not close to 0: " 
                      << g2(i) 
                      << std::endl;
        }
        if (abs(g3(i)) >= 1e-4) {
            std::cout << "        DigitCustomizedConstraints.cpp: swing foot pitch at time instance " 
                      << i 
                      << " is not close to 0: " 
                      << g3(i) 
                      << std::endl;
        }
        if (abs(g4(i)) >= 1e-4) {
            std::cout << "        DigitCustomizedConstraints.cpp: swing foot yaw at time instance " 
                      << i 
                      << " is not close to 0: " 
                      << g4(i) 
                      << std::endl;
        }
    }

    // (4) torso height always larger than 1 meter
    //     roll and pitch always close to 0
    //     yaw always close to 0 when walking forward
    for (int i = 0; i < trajPtr_->N; i++) {
        if (g6(i) <= g6_lb(i) - 1e-4) {
            std::cout << "        DigitCustomizedConstraints.cpp: torso height at time instance " 
                      << i 
                      << " is below lower limit: " 
                      << g6(i) 
                      << " < " 
                      << g6_lb(i) 
                      << std::endl;
        }
        if (abs(g9(i)) >= gp.eps_torso_angle) {
            std::cout << "        DigitCustomizedConstraints.cpp: torso roll at time instance " 
                      << i 
                      << " is too large: " 
                      << g9(i) 
                      << std::endl;
        }
        if (abs(g10(i)) >= gp.eps_torso_angle) {
            std::cout << "        DigitCustomizedConstraints.cpp: torso pitch at time instance " 
                      << i 
                      << " is too large: " 
                      << g10(i) 
                      << std::endl;
        }
        if (abs(g11(i)) >= gp.eps_torso_angle) {
            std::cout << "        DigitCustomizedConstraints.cpp: torso yaw at time instance " 
                      << i 
                      << " is too large: " 
                      << g11(i) 
                      << std::endl;
        }
    }
}

}; // namespace Digit
}; // namespace RAPTOR