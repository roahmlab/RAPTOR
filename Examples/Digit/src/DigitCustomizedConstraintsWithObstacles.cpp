#include "DigitCustomizedConstraintsWithObstacles.h"

namespace IDTO {
namespace Digit {

DigitCustomizedConstraintsWithObstacles::DigitCustomizedConstraintsWithObstacles(const Model& model_input,
                                                                                 const Eigen::VectorXi& jtype_input,
                                                                                 std::shared_ptr<Trajectories>& trajPtr_input,
                                                                                 std::shared_ptr<DynamicsConstraints>& dcPtr_input,
                                                                                 const Eigen::Array<Vec3, 1, Eigen::Dynamic>& zonotopeCenters_input,
                                                                                 const Eigen::Array<MatX, 1, Eigen::Dynamic>& zonotopeGenerators_input,
                                                                                 const VecX& q_act0_input,
                                                                                 const VecX& q_act_d0_input) : 
    jtype(jtype_input),
    trajPtr_(trajPtr_input),
    dcPtr_(dcPtr_input),
    q_act0(q_act0_input),
    q_act_d0(q_act_d0_input){
    modelPtr_ = std::make_unique<Model>(model_input);
    fkhofPtr_ = std::make_unique<ForwardKinematicsHighOrderDerivative>();
    collisionAvoidancePtr_ = std::make_shared<ZonotopeCollisionAvoidance>(zonotopeCenters_input, 
                                                                          zonotopeGenerators_input);

    if (modelPtr_->existJointName("right_toe_roll")) {
        swingfoot_id = modelPtr_->getJointId("right_toe_roll");
    }
    else {
        throw std::runtime_error("Can not find joint: right_toe_roll");
    }

    // This is right foot end transform
    // This only applies when stance foot is left foot!!!
    swingfoot_endT.R << 0,             -1, 0,
                        0.5,           0,  -sin(M_PI / 3),
                        sin(M_PI / 3), 0,  0.5;
    swingfoot_endT.p << 0, 0.05456, -0.0315;

    jointTJ = MatX::Zero(6, modelPtr_->nv);
    q = MatX::Zero(modelPtr_->nv, trajPtr_->N);
    swingfoot_xyzrpy = MatX::Zero(6, trajPtr_->N);
    pq_pz.resize(1, trajPtr_->N);
    pswingfoot_xyzrpy_pz.resize(1, trajPtr_->N);

    m = trajPtr_->N * 6 + trajPtr_->N * collisionAvoidancePtr_->numObstacles + trajPtr_->Nact * 2;

    g = VecX::Zero(m);
    g_lb = VecX::Zero(m);
    g_ub = VecX::Zero(m);
    pg_pz.resize(m, trajPtr_->varLength);
}

void DigitCustomizedConstraintsWithObstacles::compute(const VecX& z, bool compute_derivatives) {
    trajPtr_->compute(z, compute_derivatives, compute_hessian);

    // compute full joint trajectories and swing foot forward kinematics
    for (int i = 0; i < trajPtr_->N; i++) {
        VecX qi(modelPtr_->nq);
        dcPtr_->fill_independent_vector(qi, trajPtr_->q(i));
        dcPtr_->setupJointPosition(qi, compute_derivatives);
        q.col(i) = qi;

        fkhofPtr_->fk(jointT, *modelPtr_, jtype, swingfoot_id, 0, qi, swingfoot_endT, startT);

        swingfoot_xyzrpy.col(i) = fkhofPtr_->Transform2xyzrpy(jointT);

        if (compute_derivatives) {
            pq_pz(i).resize(modelPtr_->nv, trajPtr_->varLength);
            // fill in independent joints derivatives directly
            for (int j = 0; j < dcPtr_->numIndependentJoints; j++) {
                int indenpendentJointIndex = dcPtr_->return_independent_joint_index(j);
                pq_pz(i).row(indenpendentJointIndex) = trajPtr_->pq_pz(i).row(j);
            }
            // compute and fill in dependent joints derivatives
            MatX pq_dep_pz = dcPtr_->pq_dep_pq_indep * trajPtr_->pq_pz(i);
            for (int j = 0; j < dcPtr_->numDependentJoints; j++) {
                int denpendentJointIndex = dcPtr_->return_dependent_joint_index(j);
                pq_pz(i).row(denpendentJointIndex) = pq_dep_pz.row(j);
            }

            fkhofPtr_->fk_jacobian(dTdq, *modelPtr_, jtype, swingfoot_id, 0, qi, swingfoot_endT, startT);
            fkhofPtr_->Transform2xyzrpyJacobian(jointTJ, jointT, dTdq);
            pswingfoot_xyzrpy_pz(i) = jointTJ * pq_pz(i);
        }
    }

    // (1) swingfoot height always larger than 0
    //               height equal to 0 at the beginning and at the end
    //               height higher than the desired value in the middle
    VecX g1 = swingfoot_xyzrpy.row(2).transpose();

    // (2) swingfoot always flat
    VecX g2 = swingfoot_xyzrpy.row(3).transpose(); // swingfoot roll
    VecX g3 = swingfoot_xyzrpy.row(4).transpose(); // swingfoot pitch

    // (3) torso height always larger than 1 meter
    //           roll and pitch always close to 0
    VecX g6 = q.row(2).transpose(); // torso height
    VecX g7 = q.row(3).transpose(); // torso roll
    VecX g8 = q.row(4).transpose(); // torso pitch

    // (4) swing foot stays away from obstacles
    const int numberObstacles = collisionAvoidancePtr_->numObstacles;
    VecX g9 = VecX::Zero(trajPtr_->N * numberObstacles);
    MatX pg9_pz = MatX::Zero(trajPtr_->N * numberObstacles, trajPtr_->varLength);
    for (int i = 0; i < trajPtr_->N; i++) {
        const Vec3& sphereCenters = swingfoot_xyzrpy.col(i).head(3);
        collisionAvoidancePtr_->computeDistance(sphereCenters);
        g9.segment(i * numberObstacles, numberObstacles) = collisionAvoidancePtr_->distances.array() - swingfoot_sphere_radius;

        // if (i == 0) {
        //     std::cout << sphereCenters.transpose() << std::endl;
        //     std::cout << collisionAvoidancePtr_->distances << std::endl;
        // }

        if (compute_derivatives) {
            const MatX& psphereCenters_pz = pswingfoot_xyzrpy_pz(i).topRows(3);
            collisionAvoidancePtr_->computeDistance(sphereCenters, psphereCenters_pz);
            pg9_pz.block(i * numberObstacles, 0, numberObstacles, trajPtr_->varLength) = collisionAvoidancePtr_->pdistances_pz;
        }
    }

    // (5) initial condition
    VecX g10 = trajPtr_->q(0);
    VecX g11 = trajPtr_->q_d(0);

    g << g1, g2, g3, g6, g7, g8, g9, g10, g11;

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
            pg_pz.row(iter++) = pq_pz(i).row(2);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pq_pz(i).row(3);
        }
        for (int i = 0; i < trajPtr_->N; i++) {
            pg_pz.row(iter++) = pq_pz(i).row(4);
        }

        pg_pz.block(trajPtr_->N * 6, 0, trajPtr_->N * numberObstacles, trajPtr_->varLength) = pg9_pz;
    }
}

void DigitCustomizedConstraintsWithObstacles::compute_bounds() {
    // (1) swingfoot height always larger than 0
    //               height equal to 0 at the beginning and at the end
    //               height higher than the desired value in the middle
    VecX g1_lb = VecX::Zero(trajPtr_->N);
    VecX g1_ub = VecX::Constant(trajPtr_->N, 1e19);
    g1_ub(0) = 0;
    g1_ub(trajPtr_->N - 1) = 0;

    // (2) swingfoot always flat
    VecX g2_lb = VecX::Zero(trajPtr_->N);
    VecX g2_ub = VecX::Zero(trajPtr_->N);
    VecX g3_lb = VecX::Zero(trajPtr_->N);
    VecX g3_ub = VecX::Zero(trajPtr_->N);

    // (3) torso height always larger than 1 meter
    //           roll and pitch always close to 0
    VecX g6_lb = VecX::Constant(trajPtr_->N, 1);
    VecX g6_ub = VecX::Constant(trajPtr_->N, 1e19);
    VecX g7_lb = VecX::Constant(trajPtr_->N, -eps_torso_angle);
    VecX g7_ub = VecX::Constant(trajPtr_->N, eps_torso_angle);
    VecX g8_lb = VecX::Constant(trajPtr_->N, -eps_torso_angle);
    VecX g8_ub = VecX::Constant(trajPtr_->N, eps_torso_angle);

    // (4) swing foot stays away from obstacles
    VecX g9_lb = VecX::Zero(trajPtr_->N * collisionAvoidancePtr_->numObstacles);
    VecX g9_ub = VecX::Constant(trajPtr_->N * collisionAvoidancePtr_->numObstacles, 1e19);

    // (5) initial condition
    VecX g10_lb = q_act0;
    VecX g10_ub = q_act0;
    VecX g11_lb = q_act_d0;
    VecX g11_ub = q_act_d0;

    g_lb << g1_lb, g2_lb, g3_lb, g6_lb, g7_lb, g8_lb, g9_lb, g10_lb, g11_lb;
    g_ub << g1_ub, g2_ub, g3_ub, g6_ub, g7_ub, g8_ub, g9_ub, g10_ub, g11_ub;
}

}; // namespace Digit
}; // namespace IDTO