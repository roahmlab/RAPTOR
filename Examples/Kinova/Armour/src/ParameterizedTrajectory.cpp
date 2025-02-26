#include "ParameterizedTrajectory.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

BezierCurveInterval::BezierCurveInterval() {
    q0 = VecX::Zero(NUM_FACTORS);
    q_d0 = VecX::Zero(NUM_FACTORS);
    q_dd0 = VecX::Zero(NUM_FACTORS);
    Tqd0 = VecX::Zero(NUM_FACTORS);
    TTqdd0 = VecX::Zero(NUM_FACTORS);

    // compute a uniform time partition
    num_time_steps = DEFAULT_NUM_TIME_STEPS;
    s_intervals.reserve(num_time_steps + 1);
    for (int i = 0; i <= num_time_steps; i++) {
        s_intervals.push_back(i / (double)num_time_steps);
    }

    // pre-allocate memory
    q_des.resize(NUM_FACTORS, num_time_steps);
    qd_des.resize(NUM_FACTORS, num_time_steps);
    qda_des.resize(NUM_FACTORS, num_time_steps);
    qdda_des.resize(NUM_FACTORS, num_time_steps);

    ds = 1.0 / num_time_steps;

    sample_eigenvalues();
}

BezierCurveInterval::BezierCurveInterval(const VecX& q0_inp, 
                                         const VecX& q_d0_inp, 
                                         const VecX& q_dd0_inp,
                                         const VecX& k_center_inp,
                                         const VecX& k_range_inp,
                                         const double duration_inp,
                                         const std::shared_ptr<RobotInfo>& robotInfoPtr_inp,
                                         const size_t num_time_steps_inp) : 
    q0(q0_inp),
    q_d0(q_d0_inp),
    q_dd0(q_dd0_inp),
    k_center(k_center_inp),
    k_range(k_range_inp),
    duration(duration_inp),
    robotInfoPtr_(robotInfoPtr_inp),
    num_time_steps(num_time_steps_inp) {
    Tqd0 = q_d0 * duration; 
    TTqdd0 = q_dd0 * duration * duration; 

    // compute a uniform time partition
    s_intervals.reserve(num_time_steps + 1);
    for (int i = 0; i <= num_time_steps; i++) {
        s_intervals.push_back(i / (double)num_time_steps);
    }

    // pre-allocate memory
    q_des.resize(NUM_FACTORS, num_time_steps);
    qd_des.resize(NUM_FACTORS, num_time_steps);
    qda_des.resize(NUM_FACTORS, num_time_steps);
    qdda_des.resize(NUM_FACTORS, num_time_steps);

    // initialize the extrema of the k independent part of q_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        q_des_k_indep_extrema_1[i] = (2*Tqd0[i] + TTqdd0[i] + std::sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 120*k_center[i]*Tqd0[i]))/(5*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        q_des_k_indep_extrema_2[i] = (2*Tqd0[i] + TTqdd0[i] - std::sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 120*k_center[i]*Tqd0[i]))/(5*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        q_des_k_indep_extremum_1[i] = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], q_des_k_indep_extrema_1[i], duration);
        q_des_k_indep_extremum_2[i] = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], q_des_k_indep_extrema_2[i], duration);
    }

    // initialize the extrema of the k independent part of qd_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        qd_des_k_indep_extrema_1[i] = (18*Tqd0[i] + 4*TTqdd0[i] + std::sqrt(6*(54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 20*TTqdd0[i]*k_center[i] - 180*Tqd0[i]*k_center[i] + 150*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qd_des_k_indep_extrema_2[i] = (18*Tqd0[i] + 4*TTqdd0[i] - std::sqrt(6*(54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 20*TTqdd0[i]*k_center[i] - 180*Tqd0[i]*k_center[i] + 150*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qd_des_k_indep_extremum_1[i] = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qd_des_k_indep_extrema_1[i], duration);
        qd_des_k_indep_extremum_2[i] = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qd_des_k_indep_extrema_2[i], duration);
    }

    // initialize the extrema of the k independent part of qdd_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        qdd_des_k_indep_extrema_1[i] = (32*Tqd0[i] + 6*TTqdd0[i] + std::sqrt(2*(152*pow(Tqd0[i],2) + 42*Tqd0[i]*TTqdd0[i] + 3*pow(TTqdd0[i],2) - 80*TTqdd0[i]*k_center[i] - 600*Tqd0[i]*k_center[i] + 600*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qdd_des_k_indep_extrema_2[i] = (32*Tqd0[i] + 6*TTqdd0[i] - std::sqrt(2*(152*pow(Tqd0[i],2) + 42*Tqd0[i]*TTqdd0[i] + 3*pow(TTqdd0[i],2) - 80*TTqdd0[i]*k_center[i] - 600*Tqd0[i]*k_center[i] + 600*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qdd_des_k_indep_extremum_1[i] = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qdd_des_k_indep_extrema_1[i], duration);
        qdd_des_k_indep_extremum_2[i] = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qdd_des_k_indep_extrema_2[i], duration);
    }

    ds = 1.0 / num_time_steps;

    sample_eigenvalues();
}

void BezierCurveInterval::sample_eigenvalues(size_t num_samples) {
    pinocchio::Model model = robotInfoPtr_->model;

    if (robotInfoPtr_->num_joints > robotInfoPtr_->num_motors) {
        pinocchio::Model model_reduced;
        std::vector<pinocchio::JointIndex> list_of_joints_to_lock_by_id = {(pinocchio::JointIndex)model.nv};
        pinocchio::buildReducedModel(model, list_of_joints_to_lock_by_id, Eigen::VectorXd::Zero(model.nv), model_reduced);
        model_reduced.armature = model.armature.head(robotInfoPtr_->num_motors);
        model = model_reduced;
    }
    pinocchio::Data data(model);

    VecX q = VecX::Zero(robotInfoPtr_->num_motors);
    VecX qd = VecX::Zero(robotInfoPtr_->num_motors);
    VecX qdd = VecX::Zero(robotInfoPtr_->num_motors);
    
    double M_max = 0.0;
    double M_min = std::numeric_limits<double>::max();

    std::cout << "Sampling " << num_samples << " configurations..." << std::endl;
    for (int h = 0; h < num_samples; h++) {
        // generate random configurations along the trajectory
        double t = std::rand() / (double)RAND_MAX * duration;
        VecX k = VecX::Random(robotInfoPtr_->num_motors);
        computeTrajectories(k, t, q, qd, qdd);
        q += Eigen::VectorXd::Random(robotInfoPtr_->num_motors) * 0.05;

        // compute the inertia mass matrix
        MatX M = pinocchio::crba(model, data, q);
        for (int i = 0; i < M.rows(); i++) {
            for (int j = i + 1; j < M.cols(); j++) {
                M(j, i) = M(i, j);
            }
        }

        // compute the eigenvalues
        Eigen::EigenSolver<MatX> es(M);
        Eigen::VectorXd eigenvalues = es.eigenvalues().real();

        // update the maximum and minimum eigenvalues
        for (int i = 0; i < eigenvalues.size(); i++) {
            M_max = std::max(M_max, eigenvalues(i));
            M_min = std::min(M_min, eigenvalues(i));
        }
    }

    std::cout << "M_max: " << M_max << std::endl;
    std::cout << "M_min: " << M_min << std::endl;

    auto& ultimate_bound_info = robotInfoPtr_->ultimate_bound_info;
    ultimate_bound_info.M_max = M_max;
    ultimate_bound_info.M_min = M_min;

    // recompute the ultimate bound
    robotInfoPtr_->ultimate_bound_info.eps = std::sqrt(2 * ultimate_bound_info.V_m / ultimate_bound_info.M_min);
    ultimate_bound_info.qe = ultimate_bound_info.eps / ultimate_bound_info.Kr;
    ultimate_bound_info.qde = 2 * ultimate_bound_info.eps;
    ultimate_bound_info.qdae = ultimate_bound_info.eps;
    ultimate_bound_info.qddae = 2 * ultimate_bound_info.Kr * ultimate_bound_info.eps;

    std::cout << "Ultimate bound information updated:" << std::endl;
    std::cout << "eps: " << ultimate_bound_info.eps << std::endl;
    std::cout << "qe: " << Utils::rad2deg(ultimate_bound_info.qe) << " degree" << std::endl;
    std::cout << "qde: " << Utils::rad2deg(ultimate_bound_info.qde) << " degree/s" << std::endl;
    std::cout << "qdae: " << Utils::rad2deg(ultimate_bound_info.qdae) << " degree/s" << std::endl;
    std::cout << "qddae: " << Utils::rad2deg(ultimate_bound_info.qddae) << " degree/s^2" << std::endl;
}

void BezierCurveInterval::setTrajectoryParameters(const VecX& q0_inp, 
                                                  const VecX& q_d0_inp, 
                                                  const VecX& q_dd0_inp,
                                                  const VecX& k_center_inp,
                                                  const VecX& k_range_inp,
                                                  const double duration_inp) {
    q0 = q0_inp;
    q_d0 = q_d0_inp;
    q_dd0 = q_dd0_inp;
    k_center = k_center_inp;
    k_range = k_range_inp;
    duration = duration_inp;

    Tqd0 = q_d0 * duration; 
    TTqdd0 = q_dd0 * duration * duration;   

    // initialize the extrema of the k independent part of q_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        q_des_k_indep_extrema_1[i] = (2*Tqd0[i] + TTqdd0[i] + std::sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 120*k_center[i]*Tqd0[i]))/(5*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        q_des_k_indep_extrema_2[i] = (2*Tqd0[i] + TTqdd0[i] - std::sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 120*k_center[i]*Tqd0[i]))/(5*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        q_des_k_indep_extremum_1[i] = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], q_des_k_indep_extrema_1[i], duration);
        q_des_k_indep_extremum_2[i] = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], q_des_k_indep_extrema_2[i], duration);
    }

    // initialize the extrema of the k independent part of qd_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        qd_des_k_indep_extrema_1[i] = (18*Tqd0[i] + 4*TTqdd0[i] + std::sqrt(6*(54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 20*TTqdd0[i]*k_center[i] - 180*Tqd0[i]*k_center[i] + 150*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qd_des_k_indep_extrema_2[i] = (18*Tqd0[i] + 4*TTqdd0[i] - std::sqrt(6*(54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 20*TTqdd0[i]*k_center[i] - 180*Tqd0[i]*k_center[i] + 150*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qd_des_k_indep_extremum_1[i] = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qd_des_k_indep_extrema_1[i], duration);
        qd_des_k_indep_extremum_2[i] = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qd_des_k_indep_extrema_2[i], duration);
    }

    // initialize the extrema of the k independent part of qdd_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        qdd_des_k_indep_extrema_1[i] = (32*Tqd0[i] + 6*TTqdd0[i] + std::sqrt(2*(152*pow(Tqd0[i],2) + 42*Tqd0[i]*TTqdd0[i] + 3*pow(TTqdd0[i],2) - 80*TTqdd0[i]*k_center[i] - 600*Tqd0[i]*k_center[i] + 600*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qdd_des_k_indep_extrema_2[i] = (32*Tqd0[i] + 6*TTqdd0[i] - std::sqrt(2*(152*pow(Tqd0[i],2) + 42*Tqd0[i]*TTqdd0[i] + 3*pow(TTqdd0[i],2) - 80*TTqdd0[i]*k_center[i] - 600*Tqd0[i]*k_center[i] + 600*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qdd_des_k_indep_extremum_1[i] = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qdd_des_k_indep_extrema_1[i], duration);
        qdd_des_k_indep_extremum_2[i] = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qdd_des_k_indep_extrema_2[i], duration);
    }                                          
}

void BezierCurveInterval::computeTrajectories(
    const VecX& k,
    const double t,
    VecX& q,
    VecX& qd,
    VecX& qdd) const {
    assert(t >= 0 && t <= duration);
    assert(k.size() == NUM_FACTORS);
    assert(q.size() >= NUM_FACTORS);
    assert(qd.size() >= NUM_FACTORS);
    assert(qdd.size() >= NUM_FACTORS);
    const double s = t / duration;
    for (int i = 0; i < NUM_FACTORS; i++) {
        const double k_actual = k_center[i] + k_range[i] * k[i];
        q[i] = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, s);
        qd[i] = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, s, duration);
        qdd[i] = qdd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, s, duration);
    }
}

void BezierCurveInterval::makePolyZono(const int s_ind) {
    assert(s_ind < num_time_steps);

    const auto& ultimate_bound_info = robotInfoPtr_->ultimate_bound_info;

    // const double s_lb = s_ind * ds;
    // const double s_ub = (s_ind + 1) * ds;
    const double s_lb = s_intervals[s_ind];
    const double s_ub = s_intervals[s_ind + 1];

    const Interval t_int(s_lb, s_ub);

    for (int i = 0; i < NUM_FACTORS; i++) {
        const double k_range_elt = k_range(i);

        // Part 1: q_des
        double k_dep_coeff_lb = pow(s_lb,3) * (6 * pow(s_lb,2) - 15 * s_lb + 10);
        double k_dep_coeff_ub = pow(s_ub,3) * (6 * pow(s_ub,2) - 15 * s_ub + 10);
        double k_dep_coeff_center = (k_dep_coeff_ub + k_dep_coeff_lb) * 0.5 * k_range_elt;
        double k_dep_coeff_radius = (k_dep_coeff_ub - k_dep_coeff_lb) * 0.5 * k_range_elt;
        
        double k_indep_lb = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], s_lb, duration);
        double k_indep_ub = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], s_ub, duration);
        if (k_indep_lb > k_indep_ub) {
            std::swap(k_indep_lb, k_indep_ub);
        }
        if (s_lb < q_des_k_indep_extrema_1[i] && q_des_k_indep_extrema_1[i] < s_ub) {
            k_indep_lb = std::min(k_indep_lb, q_des_k_indep_extremum_1[i]);
            k_indep_ub = std::max(k_indep_ub, q_des_k_indep_extremum_1[i]);
        }
        if (s_lb < q_des_k_indep_extrema_2[i] && q_des_k_indep_extrema_2[i] < s_ub) {
            k_indep_lb = std::min(k_indep_lb, q_des_k_indep_extremum_2[i]);
            k_indep_ub = std::max(k_indep_ub, q_des_k_indep_extremum_2[i]);
        }
        double k_indep_radius = (k_indep_ub - k_indep_lb) * 0.5;
        double q_des_center = (k_indep_lb + k_indep_ub) * 0.5;
        
        // PZSparse of q_des
        uint32_t q_des_degree[2][NUM_VARIABLES] = {0};
        double q_des_coeff[2] = {0};
        q_des_degree[0][i] = 1; // k
        q_des_coeff[0] = k_dep_coeff_center; // k
        q_des_degree[1][i + NUM_FACTORS * 1] = 1; // qe
        q_des_coeff[1] = k_dep_coeff_radius + k_indep_radius + ultimate_bound_info.qe; // qe
        PZSparse q_des_range(0, q_des_coeff, q_des_degree, 2);
        q_des(i, s_ind) = PZSparse(q_des_center, q_des_coeff, q_des_degree, 2);

        // Part 2: qd_des
        if (s_ub <= 0.5) {
            k_dep_coeff_lb = (30 * pow(s_lb,2) * pow(s_lb - 1,2)) / duration;
            k_dep_coeff_ub = (30 * pow(s_ub,2) * pow(s_ub - 1,2)) / duration;
        }
        else if (s_lb >= 0.5) {
            k_dep_coeff_lb = (30 * pow(s_ub,2) * pow(s_ub - 1,2)) / duration;
            k_dep_coeff_ub = (30 * pow(s_lb,2) * pow(s_lb - 1,2)) / duration;
        }
        else {
            k_dep_coeff_lb = std::min((30 * pow(s_lb,2) * pow(s_lb - 1,2)) / duration,
                                 (30 * pow(s_ub,2) * pow(s_ub - 1,2)) / duration);
            k_dep_coeff_ub = (30 * pow(0.5,2) * pow(0.5 - 1,2)) / duration;
        }

        k_dep_coeff_center = (k_dep_coeff_ub + k_dep_coeff_lb) * 0.5 * k_range_elt; // Have to scale to [-1,1] in order to fit in PZSparse
        k_dep_coeff_radius = (k_dep_coeff_ub - k_dep_coeff_lb) * 0.5 * k_range_elt; // Have to scale to [-1,1] in order to fit in PZSparse

        k_indep_lb = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], s_lb, duration);
        k_indep_ub = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], s_ub, duration);
        if (k_indep_lb > k_indep_ub) {
            std::swap(k_indep_lb, k_indep_ub);
        }
        if (s_lb < qd_des_k_indep_extrema_1[i] && qd_des_k_indep_extrema_1[i] < s_ub) {
            k_indep_lb = std::min(k_indep_lb, qd_des_k_indep_extremum_1[i]);
            k_indep_ub = std::max(k_indep_ub, qd_des_k_indep_extremum_1[i]);
        }
        if (s_lb < qd_des_k_indep_extrema_2[i] && qd_des_k_indep_extrema_2[i] < s_ub) {
            k_indep_lb = std::min(k_indep_lb, qd_des_k_indep_extremum_2[i]);
            k_indep_ub = std::max(k_indep_ub, qd_des_k_indep_extremum_2[i]);
        }
        k_indep_radius = (k_indep_ub - k_indep_lb) * 0.5;
        double qd_des_center = (k_indep_lb + k_indep_ub) * 0.5;

        double qd_des_coeff[] = {k_dep_coeff_center, k_dep_coeff_radius + k_indep_radius + ultimate_bound_info.qde};

        uint32_t qd_des_degree[2][NUM_VARIABLES] = {0};
        qd_des_degree[0][i] = 1; // k
        qd_des_degree[1][i + NUM_FACTORS * 2] = 1; // qde

        // qd_des_int = qd_des_center + qd_des_coeff[0] * k + qd_des_coeff[1] * qde;
        qd_des(i, s_ind) = PZSparse(qd_des_center, qd_des_coeff, qd_des_degree, 2);

        double qda_des_coeff[] = {k_dep_coeff_center, k_dep_coeff_radius + k_indep_radius + ultimate_bound_info.qdae};

        uint32_t qda_des_degree[2][NUM_VARIABLES] = {0};
        qda_des_degree[0][i] = 1; // k
        qda_des_degree[1][i + NUM_FACTORS * 3] = 1; // qdae

        // qda_des_int = qd_des_center + qda_des_coeff[0] * k + qda_des_coeff[1] * qdae;
        qda_des(i, s_ind) = PZSparse(qd_des_center, qda_des_coeff, qda_des_degree, 2);

        // Part 3: qdd_des
        double temp_lb = (60 * s_lb * (2 * pow(s_lb,2) - 3 * s_lb + 1)) / duration / duration;
        double temp_ub = (60 * s_ub * (2 * pow(s_ub,2) - 3 * s_ub + 1)) / duration / duration;
        if (s_ub <= QDD_DES_K_DEP_MAXIMA) { // monotonically increasing region
            k_dep_coeff_lb = temp_lb;
            k_dep_coeff_ub = temp_ub;
        }
        else if (s_lb <= QDD_DES_K_DEP_MAXIMA) { // maxima lives inside
            k_dep_coeff_lb = std::min(temp_lb, temp_ub);
            k_dep_coeff_ub = (60 * QDD_DES_K_DEP_MAXIMA * (2 * pow(QDD_DES_K_DEP_MAXIMA,2) - 3 * QDD_DES_K_DEP_MAXIMA + 1)) / duration / duration;
        }
        else if (s_ub <= QDD_DES_K_DEP_MINIMA) { // monotonically decreasing region
            k_dep_coeff_lb = temp_ub;   
            k_dep_coeff_ub = temp_lb;
        }
        else if (s_lb <= QDD_DES_K_DEP_MINIMA) { // minima lives inside
            k_dep_coeff_lb = (60 * QDD_DES_K_DEP_MINIMA * (2 * pow(QDD_DES_K_DEP_MINIMA,2) - 3 * QDD_DES_K_DEP_MINIMA + 1)) / duration / duration;
            k_dep_coeff_ub = std::max(temp_lb, temp_ub);
        }
        else { // monotonically increasing region
            k_dep_coeff_lb = temp_lb;
            k_dep_coeff_ub = temp_ub;
        }
        
        k_dep_coeff_center = (k_dep_coeff_ub + k_dep_coeff_lb) * 0.5 * k_range_elt; // Have to scale to [-1,1] in order to fit in PZSparse
        k_dep_coeff_radius = (k_dep_coeff_ub - k_dep_coeff_lb) * 0.5 * k_range_elt;

        k_indep_lb = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], s_lb, duration);
        k_indep_ub = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], s_ub, duration);
        if (k_indep_lb > k_indep_ub) {
            std::swap(k_indep_lb, k_indep_ub);
        }
        if (s_lb < qdd_des_k_indep_extrema_1[i] && qdd_des_k_indep_extrema_1[i] < s_ub) {
            k_indep_lb = std::min(k_indep_lb, qdd_des_k_indep_extremum_1[i]);
            k_indep_ub = std::max(k_indep_ub, qdd_des_k_indep_extremum_1[i]);
        }
        if (s_lb < qdd_des_k_indep_extrema_2[i] && qdd_des_k_indep_extrema_2[i] < s_ub) {
            k_indep_lb = std::min(k_indep_lb, qdd_des_k_indep_extremum_2[i]);
            k_indep_ub = std::max(k_indep_ub, qdd_des_k_indep_extremum_2[i]);
        }
        k_indep_radius = (k_indep_ub - k_indep_lb) * 0.5;
        double qdd_des_center = (k_indep_lb + k_indep_ub) * 0.5;

        double qdd_des_coeff[] = {k_dep_coeff_center, k_dep_coeff_radius + k_indep_radius + ultimate_bound_info.qddae};

        uint32_t qdd_des_degree[2][NUM_VARIABLES] = {0};
        qdd_des_degree[0][i] = 1; // k
        qdd_des_degree[1][i + NUM_FACTORS * 4] = 1; // qddae

        // qdd_des_int = qdd_des_center + qdd_des_coeff[0] * k + qdd_des_coeff[1] * qdde;
        qdda_des(i, s_ind) = PZSparse(qdd_des_center, qdd_des_coeff, qdd_des_degree, 2);
    }
}

void BezierCurveInterval::returnJointPositionExtremum(double* extremum, const double* k) const {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        double k_actual = k_center(i) + k_range(i) * k[i];

        // list all possible extremas
        double extrema1 = 0;
        double extrema2 = (2*Tqd0[i] + TTqdd0[i] + std::sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] - 120*k_actual*Tqd0[i] + pow(TTqdd0[i],2))) / (5*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema3 = (2*Tqd0[i] + TTqdd0[i] - std::sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] - 120*k_actual*Tqd0[i] + pow(TTqdd0[i],2))) / (5*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema4 = 1;

        // get extremums of all extremas
        double extremum1 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema1);
        double extremum2 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema2);
        double extremum3 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema3);
        double extremum4 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema4);

        // find the min and max values
        double minPosition = std::min(extremum1, extremum4);
        double maxPosition = std::max(extremum1, extremum4);
        if (0 <= extrema2 && extrema2 <= 1) { // check if this extrema is inside the time range [0,1]
            minPosition = std::min(minPosition, extremum2);
            maxPosition = std::max(maxPosition, extremum2);
        }
        if (0 <= extrema3 && extrema3 <= 1) { // check if this extrema is inside the time range [0,1]
            minPosition = std::min(minPosition, extremum3);
            maxPosition = std::max(maxPosition, extremum3);
        }

        extremum[i              ] = minPosition;
        extremum[i + NUM_FACTORS] = maxPosition;
    }
}

void BezierCurveInterval::returnJointPositionExtremumGradient(double* extremumGradient, const double* k) const {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        double k_actual = k_center(i) + k_range(i) * k[i];

        // list all possible extremas
        double extrema1 = 0;
        double extrema2 = (2*Tqd0[i] + TTqdd0[i] + std::sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] - 120*k_actual*Tqd0[i] + pow(TTqdd0[i],2))) / (5*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema3 = (2*Tqd0[i] + TTqdd0[i] - std::sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] - 120*k_actual*Tqd0[i] + pow(TTqdd0[i],2))) / (5*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema4 = 1;

        // get extremums of all extremas
        double extremum1 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema1);
        double extremum2 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema2);
        double extremum3 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema3);
        double extremum4 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema4);

        // find the min and max values
        double minPosition;
        int minId;
        double maxPosition;
        int maxId;

        if (extremum1 < extremum4) {
            minPosition = extremum1;
            minId = 1;

            maxPosition = extremum4;
            maxId = 4;
        }
        else {
            minPosition = extremum4;
            minId = 4;

            maxPosition = extremum1;
            maxId = 1;
        }

        if (0 <= extrema2 && extrema2 <= 1) { // check if this extrema is inside the time range [0,1]
            if (extremum2 < minPosition) {
                minPosition = extremum2;
                minId = 2;
            }
            if (maxPosition < extremum2) {
                maxPosition = extremum2;
                maxId = 2;
            }
        }
        if (0 <= extrema3 && extrema3 <= 1) { // check if this extrema is inside the time range [0,1]
            if (extremum3 < minPosition) {
                minPosition = extremum3;
                minId = 3;
            }
            if (maxPosition < extremum3) {
                maxPosition = extremum3;
                maxId = 3;
            }
        }

        double minPositionGradient;
        double maxPositionGradient;

        switch (minId) {
            case 1: // t = 0
                minPositionGradient = 0.0;
                break;
            case 2: // t = extrema2
                minPositionGradient = q_des_extrema2_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 3: // t = extrema3
                minPositionGradient = q_des_extrema3_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 4: // t = 1
                minPositionGradient = 1.0;
                break;
            default:
                break;
        }

        switch (maxId) {
            case 1: // t = 0
                maxPositionGradient = 0.0;
                break;
            case 2: // t = extrema2
                maxPositionGradient = q_des_extrema2_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 3: // t = extrema3
                maxPositionGradient = q_des_extrema3_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 4: // t = 1
                maxPositionGradient = 1.0;
                break;
            default:
                break;
        }

        for (int j = 0; j < NUM_FACTORS; j++) {
            if (i == j) {
                extremumGradient[(i              ) * NUM_FACTORS + j] = minPositionGradient * k_range(i);
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = maxPositionGradient * k_range(i);
            }
            else {
                extremumGradient[(i              ) * NUM_FACTORS + j] = 0.0;
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = 0.0;
            }
        }
    }
}

void BezierCurveInterval::returnJointVelocityExtremum(double* extremum, const double* k) const {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        double k_actual = k_center(i) + k_range(i) * k[i];

        // list all possible extremas
        double extrema1 = 0;
        double extrema2 = (18*Tqd0[i] - 30*k_actual + 4*TTqdd0[i] + std::sqrt(6*(150*pow(k_actual,2) - 180*k_actual*Tqd0[i] - 20*k_actual*TTqdd0[i] + 54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema3 = (18*Tqd0[i] - 30*k_actual + 4*TTqdd0[i] - std::sqrt(6*(150*pow(k_actual,2) - 180*k_actual*Tqd0[i] - 20*k_actual*TTqdd0[i] + 54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema4 = 1;

        // get extremums of all extremas
        double extremum1 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema1, duration);
        double extremum2 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema2, duration);
        double extremum3 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema3, duration);
        double extremum4 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema4, duration);

        // find the min and max values
        double minVelocity = std::min(extremum1, extremum4);
        double maxVelocity = std::max(extremum1, extremum4);
        if (0 <= extrema2 && extrema2 <= 1) { // check if this extrema is inside the time range [0,1]
            minVelocity = std::min(minVelocity, extremum2); 
            maxVelocity = std::max(maxVelocity, extremum2);
        }
        if (0 <= extrema3 && extrema3 <= 1) { // check if this extrema is inside the time range [0,1]
            minVelocity = std::min(minVelocity, extremum3);
            maxVelocity = std::max(maxVelocity, extremum3);
        }

        extremum[i              ] = minVelocity;
        extremum[i + NUM_FACTORS] = maxVelocity;
    }
}

void BezierCurveInterval::returnJointVelocityExtremumGradient(double* extremumGradient, const double* k) const {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        double k_actual = k_center(i) + k_range(i) * k[i];

        // list all possible extremas
        double extrema1 = 0;
        double extrema2 = (18*Tqd0[i] - 30*k_actual + 4*TTqdd0[i] + std::sqrt(6*(150*pow(k_actual,2) - 180*k_actual*Tqd0[i] - 20*k_actual*TTqdd0[i] + 54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema3 = (18*Tqd0[i] - 30*k_actual + 4*TTqdd0[i] - std::sqrt(6*(150*pow(k_actual,2) - 180*k_actual*Tqd0[i] - 20*k_actual*TTqdd0[i] + 54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema4 = 1;

        // get extremums of all extremas
        double extremum1 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema1, duration);
        double extremum2 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema2, duration);
        double extremum3 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema3, duration);
        double extremum4 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema4, duration);

        // find the min and max values
        double minVelocity;
        int minId;
        double maxVelocity;
        int maxId;

        if (extremum1 < extremum4) {
            minVelocity = extremum1;
            minId = 1;

            maxVelocity = extremum4;
            maxId = 4;
        }
        else {
            minVelocity = extremum4;
            minId = 4;

            maxVelocity = extremum1;
            maxId = 1;
        }

        if (0 <= extrema2 && extrema2 <= 1) { // check if this extrema is inside the time range [0,1]
            if (extremum2 < minVelocity) {
                minVelocity = extremum2;
                minId = 2;
            }
            if (maxVelocity < extremum2) {
                maxVelocity = extremum2;
                maxId = 2;
            }
        }
        if (0 <= extrema3 && extrema3 <= 1) { // check if this extrema is inside the time range [0,1]
            if (extremum3 < minVelocity) {
                minVelocity = extremum3;
                minId = 3;
            }
            if (maxVelocity < extremum3) {
                maxVelocity = extremum3;
                maxId = 3;
            }
        }

        double minVelocityGradient;
        double maxVelocityGradient;

        switch (minId) {
            case 1: // t = 0
                minVelocityGradient = 0.0;
                break;
            case 2: // t = extrema2
                minVelocityGradient = qd_des_extrema2_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 3: // t = extrema3
                minVelocityGradient = qd_des_extrema3_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 4: // t = 1
                minVelocityGradient = 0.0;
                break;
            default:
                break;
        }

        switch (maxId) {
            case 1: // t = 0
                maxVelocityGradient = 0.0;
                break;
            case 2: // t = extrema2
                maxVelocityGradient = qd_des_extrema2_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 3: // t = extrema3
                maxVelocityGradient = qd_des_extrema3_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 4: // t = 1
                maxVelocityGradient = 0.0;
                break;
            default:
                break;
        }

        for (int j = 0; j < NUM_FACTORS; j++) {
            if (i == j) {
                extremumGradient[(i              ) * NUM_FACTORS + j] = minVelocityGradient * k_range(i) / duration;
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = maxVelocityGradient * k_range(i) / duration;
            }
            else {
                extremumGradient[(i              ) * NUM_FACTORS + j] = 0.0;
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = 0.0;
            }
        }
    }
}

double q_des_func(const double q0, const double Tqd0, const double TTqdd0, const double k, const double s) {
    double B0 = -pow(s - 1,5);
    double B1 = 5*s*pow(s - 1,4);
    double B2 = -10*pow(s,2)*pow(s - 1,3);
    double B3 = 10*pow(s,3)*pow(s - 1,2);
    double B4 = -5*pow(s,4)*(s - 1);
    double B5 = pow(s,5);
    double beta0 = q0;
    double beta1 = q0 + Tqd0/5;
    double beta2 = q0 + (2*Tqd0)/5 + TTqdd0/20;
    double beta3 = q0 + k;
    double beta4 = q0 + k;
    double beta5 = q0 + k;
    return B0 * beta0 + B1 * beta1 + B2 * beta2 + B3 * beta3 + B4 * beta4 + B5 * beta5;
}

double qd_des_func(const double q0, const double Tqd0, const double TTqdd0, const double k, const double s, const double T) {
    double dB0 = pow(s-1.0,4.0)*-5.0;
    double dB1 = s*pow(s-1.0,3.0)*2.0E+1+pow(s-1.0,4.0)*5.0;
    double dB2 = s*pow(s-1.0,3.0)*-2.0E+1-(s*s)*pow(s-1.0,2.0)*3.0E+1;
    double dB3 = pow(s,3.0)*(s*2.0-2.0)*1.0E+1+(s*s)*pow(s-1.0,2.0)*3.0E+1;
    double dB4 = pow(s,3.0)*(s-1.0)*-2.0E+1-pow(s,4.0)*5.0;
    double dB5 = pow(s,4.0)*5.0;
    double beta0 = q0;
    double beta1 = q0 + Tqd0/5;
    double beta2 = q0 + (2*Tqd0)/5 + TTqdd0/20;
    double beta3 = q0 + k;
    double beta4 = q0 + k;
    double beta5 = q0 + k;
    double traj = dB0 * beta0 + dB1 * beta1 + dB2 * beta2 + dB3 * beta3 + dB4 * beta4 + dB5 * beta5;
    return traj / T;
}

double qdd_des_func(const double q0, const double Tqd0, const double TTqdd0, const double k, const double s, const double T) {
    double t2 = s*2.0;
    double t3 = s*s;
    double t4 = s*s*s;
    double t5 = s-1.0;
    double t6 = t2-2.0;
    double t7 = t4*2.0E+1;
    double t8 = t5*t5;
    double t9 = t5*t5*t5;
    double t10 = t9*2.0E+1;
    double t11 = s*t8*6.0E+1;
    double t12 = -t10;
    double ddB0 = t12;
    double ddB1 = t9*4.0E+1+t11;
    double ddB2 = t12-s*t8*1.2E+2-t3*t6*3.0E+1;
    double ddB3 = t7+t11+t3*t6*6.0E+1;
    double ddB4 = t4*-4.0E+1-t3*t5*6.0E+1;
    double ddB5 = t7;
    double beta0 = q0;
    double beta1 = q0 + Tqd0/5;
    double beta2 = q0 + (2*Tqd0)/5 + TTqdd0/20;
    double beta3 = q0 + k;
    double beta4 = q0 + k;
    double beta5 = q0 + k;
    double traj = ddB0 * beta0 + ddB1 * beta1 + ddB2 * beta2 + ddB3 * beta3 + ddB4 * beta4 + ddB5 * beta5;
    return traj / (T * T);
}

double q_des_extrema2_k_derivative(double q0, double Tqd0, double TTqdd0, double k) {
    double t2 = k+q0;
    double t3 = Tqd0*2.0;
    double t4 = Tqd0*6.0;
    double t5 = Tqd0*Tqd0;
    double t6 = TTqdd0*TTqdd0;
    double t7 = k*1.2E+1;
    double t8 = Tqd0*TTqdd0*1.4E+1;
    double t10 = Tqd0/5.0;
    double t11 = Tqd0*(2.0/5.0);
    double t12 = k*Tqd0*1.2E+2;
    double t13 = TTqdd0/2.0E+1;
    double t9 = -t7;
    double t14 = t5*6.4E+1;
    double t15 = -t12;
    double t16 = q0+t10;
    double t18 = q0+t11+t13;
    double t17 = TTqdd0+t4+t9;
    double t24 = t6+t8+t14+t15;
    double t19 = 1.0/t17;
    double t25 = std::sqrt(t24);
    double t20 = t19*t19;
    double t21 = t19*t19*t19;
    double t23 = t19*t19*t19*t19*t19;
    double t26 = 1.0/t25;
    double t27 = TTqdd0+t3+t25;
    double t22 = t20*t20;
    double t28 = t27*t27;
    double t29 = t27*t27*t27;
    double t31 = t27*t27*t27*t27*t27;
    double t32 = Tqd0*t19*t26*1.2E+1;
    double t34 = (t19*t27)/5.0;
    double t35 = t20*t27*(1.2E+1/5.0);
    double t30 = t28*t28;
    double t33 = -t32;
    double t36 = t34-1.0;
    double t37 = t36*t36;
    double t38 = t36*t36*t36;
    double t40 = t33+t35;
    double t39 = t37*t37;
    return (t23*t31)/3.125E+3+t2*(t20*t20*t20)*t31*(1.2E+1/6.25E+2)+t21*t29*t37*(2.0/2.5E+1)-(t22*t30*t36)/1.25E+2+q0*t39*(t32-t35)*5.0+t2*t22*t29*t37*(7.2E+1/2.5E+1)-t2*t23*t30*t36*(4.8E+1/1.25E+2)+t16*t20*t27*t39*1.2E+1-t18*t21*t28*t38*(4.8E+1/5.0)+(t2*t22*t30*(t32-t35))/1.25E+2-Tqd0*t2*t23*t26*t30*(1.2E+1/1.25E+2)-Tqd0*t16*t19*t26*t39*6.0E+1-t2*t21*t29*t36*(t32-t35)*(4.0/2.5E+1)-t16*t19*t27*t38*(t32-t35)*4.0+t18*t20*t28*t37*(t32-t35)*(6.0/5.0)-Tqd0*t2*t21*t26*t28*t37*(7.2E+1/5.0)+Tqd0*t2*t22*t26*t29*t36*(4.8E+1/2.5E+1)+Tqd0*t18*t20*t26*t27*t38*4.8E+1;
}

double q_des_extrema3_k_derivative(double q0, double Tqd0, double TTqdd0, double k) {
    double t2 = k+q0;
    double t3 = Tqd0*2.0;
    double t4 = Tqd0*6.0;
    double t5 = Tqd0*Tqd0;
    double t6 = TTqdd0*TTqdd0;
    double t7 = k*1.2E+1;
    double t8 = Tqd0*TTqdd0*1.4E+1;
    double t10 = Tqd0/5.0;
    double t11 = Tqd0*(2.0/5.0);
    double t12 = k*Tqd0*1.2E+2;
    double t13 = TTqdd0/2.0E+1;
    double t9 = -t7;
    double t14 = t5*6.4E+1;
    double t15 = -t12;
    double t16 = q0+t10;
    double t18 = q0+t11+t13;
    double t17 = TTqdd0+t4+t9;
    double t24 = t6+t8+t14+t15;
    double t19 = 1.0/t17;
    double t25 = std::sqrt(t24);
    double t20 = t19*t19;
    double t21 = t19*t19*t19;
    double t23 = t19*t19*t19*t19*t19;
    double t26 = 1.0/t25;
    double t27 = -t25;
    double t22 = t20*t20;
    double t28 = TTqdd0+t3+t27;
    double t33 = Tqd0*t19*t26*1.2E+1;
    double t29 = t28*t28;
    double t30 = t28*t28*t28;
    double t32 = t28*t28*t28*t28*t28;
    double t34 = (t19*t28)/5.0;
    double t35 = t20*t28*(1.2E+1/5.0);
    double t31 = t29*t29;
    double t36 = t34-1.0;
    double t40 = t33+t35;
    double t37 = t36*t36;
    double t38 = t36*t36*t36;
    double t39 = t37*t37;
    return (t23*t32)/3.125E+3+t2*(t20*t20*t20)*t32*(1.2E+1/6.25E+2)-q0*t39*t40*5.0+t21*t30*t37*(2.0/2.5E+1)-(t22*t31*t36)/1.25E+2+t2*t22*t30*t37*(7.2E+1/2.5E+1)-t2*t23*t31*t36*(4.8E+1/1.25E+2)-(t2*t22*t31*t40)/1.25E+2+t16*t20*t28*t39*1.2E+1-t18*t21*t29*t38*(4.8E+1/5.0)+Tqd0*t2*t23*t26*t31*(1.2E+1/1.25E+2)+Tqd0*t16*t19*t26*t39*6.0E+1+t2*t21*t30*t36*t40*(4.0/2.5E+1)+t16*t19*t28*t38*t40*4.0-t18*t20*t29*t37*t40*(6.0/5.0)+Tqd0*t2*t21*t26*t29*t37*(7.2E+1/5.0)-Tqd0*t2*t22*t26*t30*t36*(4.8E+1/2.5E+1)-Tqd0*t18*t20*t26*t28*t38*4.8E+1;
}

double qd_des_extrema2_k_derivative(double q0, double Tqd0, double TTqdd0, double k) {
    double t2 = k+q0;
    double t3 = k*k;
    double t4 = Tqd0*6.0;
    double t5 = TTqdd0*4.0;
    double t6 = Tqd0*Tqd0;
    double t7 = TTqdd0*TTqdd0;
    double t8 = k*1.2E+1;
    double t9 = k*3.0E+1;
    double t10 = Tqd0*1.8E+1;
    double t11 = TTqdd0*2.0E+1;
    double t13 = Tqd0*TTqdd0*1.4E+1;
    double t14 = std::sqrt(6.0);
    double t17 = k*3.0E+2;
    double t18 = Tqd0*1.8E+2;
    double t21 = k*TTqdd0*-2.0E+1;
    double t24 = k*Tqd0*-1.8E+2;
    double t12 = k*t11;
    double t15 = -t8;
    double t16 = -t9;
    double t19 = t6*5.4E+1;
    double t20 = k*t18;
    double t22 = -t17;
    double t23 = t3*1.5E+2;
    double t25 = TTqdd0+t4+t15;
    double t31 = t11+t18+t22;
    double t32 = t7+t13+t19+t21+t23+t24;
    double t26 = 1.0/t25;
    double t33 = std::sqrt(t32);
    double t27 = t26*t26;
    double t28 = t26*t26*t26;
    double t30 = t26*t26*t26*t26*t26;
    double t34 = 1.0/t33;
    double t35 = t14*t33;
    double t29 = t27*t27;
    double t36 = t5+t10+t16+t35;
    double t40 = (t14*t31*t34)/2.0;
    double t37 = t36*t36;
    double t38 = t36*t36*t36;
    double t41 = t40+3.0E+1;
    double t42 = (t26*t36)/5.0;
    double t43 = t27*t36*(6.0/5.0);
    double t44 = (t26*t36)/1.0E+1;
    double t39 = t37*t37;
    double t45 = t42-2.0;
    double t46 = t44-1.0;
    double t49 = (t26*t41)/1.0E+1;
    double t47 = t46*t46;
    double t48 = t46*t46*t46;
    double t50 = -t49;
    double t51 = t27*t36*t48*2.4E+1;
    double t52 = t28*t37*t47*(3.6E+1/5.0);
    double t53 = t43+t50;
    double t54 = t26*t41*t48*2.0;
    double t56 = t27*t36*t41*t47*(3.0/5.0);
    double t55 = -t54;
    double t57 = -t56;
    double t58 = t26*t36*t47*t53*6.0;
    double t59 = t27*t37*t46*t53*(3.0/5.0);
    return (q0+Tqd0/5.0)*(t51+t55+t58+t48*t53*2.0E+1)+t2*(t52+t57+t59+(t28*t38*(t27*t36*(1.2E+1/5.0)-(t26*t41)/5.0))/1.0E+2+t29*t38*t45*(9.0/2.5E+1)-t28*t37*t41*t45*(3.0/1.0E+2))-t2*(t30*t39*(3.0/1.25E+2)-(t29*t38*t41)/5.0E+2+t29*t38*t46*(1.8E+1/2.5E+1)+(t28*t38*t53)/5.0E+1-t28*t37*t41*t46*(3.0/5.0E+1))-(q0+Tqd0*(2.0/5.0)+TTqdd0/2.0E+1)*(t51+t52+t55+t57+t58+t59)-q0*t48*t53*2.0E+1+t2*t30*t39*(3.0/1.25E+2)+t27*t37*t47*(3.0/1.0E+1)+(t28*t38*t45)/1.0E+2-(t28*t38*t46)/5.0E+1-(t2*t29*t38*t41)/5.0E+2;
}

double qd_des_extrema3_k_derivative(double q0, double Tqd0, double TTqdd0, double k) {
    double t2 = k+q0;
    double t3 = k*k;
    double t4 = Tqd0*6.0;
    double t5 = TTqdd0*4.0;
    double t6 = Tqd0*Tqd0;
    double t7 = TTqdd0*TTqdd0;
    double t8 = k*1.2E+1;
    double t9 = k*3.0E+1;
    double t10 = Tqd0*1.8E+1;
    double t12 = TTqdd0*2.0E+1;
    double t14 = Tqd0*TTqdd0*1.4E+1;
    double t15 = std::sqrt(6.0);
    double t17 = k*3.0E+2;
    double t19 = Tqd0*1.8E+2;
    double t22 = k*TTqdd0*-2.0E+1;
    double t25 = k*Tqd0*-1.8E+2;
    double t11 = -t5;
    double t13 = k*t12;
    double t16 = -t8;
    double t18 = -t10;
    double t20 = t6*5.4E+1;
    double t21 = k*t19;
    double t23 = -t17;
    double t24 = t3*1.5E+2;
    double t26 = TTqdd0+t4+t16;
    double t32 = t12+t19+t23;
    double t33 = t7+t14+t20+t22+t24+t25;
    double t27 = 1.0/t26;
    double t34 = std::sqrt(t33);
    double t28 = t27*t27;
    double t29 = t27*t27*t27;
    double t31 = t27*t27*t27*t27*t27;
    double t35 = 1.0/t34;
    double t36 = t15*t34;
    double t30 = t28*t28;
    double t37 = t9+t11+t18+t36;
    double t38 = pow(t5-t9+t10-t36,2.0);
    double t39 = -pow(t5-t9+t10-t36,3.0);
    double t41 = (t15*t32*t35)/2.0;
    double t43 = t27*(t5-t9+t10-t36)*(-1.0/5.0);
    double t44 = t28*(t5-t9+t10-t36)*(-6.0/5.0);
    double t45 = t27*(t5-t9+t10-t36)*(-1.0/1.0E+1);
    double t48 = pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,2.0);
    double t49 = -pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0);
    double t52 = t28*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*(t5-t9+t10-t36)*2.4E+1;
    double t40 = t38*t38;
    double t42 = t41-3.0E+1;
    double t46 = t43+2.0;
    double t47 = t45+1.0;
    double t53 = t29*t38*t48*(3.6E+1/5.0);
    double t50 = (t27*t42)/1.0E+1;
    double t55 = t27*t42*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*-2.0;
    double t56 = t27*t42*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*2.0;
    double t57 = t28*t42*t48*(t5-t9+t10-t36)*(-3.0/5.0);
    double t58 = t28*t42*t48*(t5-t9+t10-t36)*(3.0/5.0);
    double t51 = -t50;
    double t59 = t27*t48*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*(t5-t9+t10-t36)*6.0;
    double t60 = t28*t38*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*(3.0/5.0);
    double t54 = t44+t51;
    return (q0+Tqd0/5.0)*(t52+t56+t59+pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*2.0E+1)+t2*(t53+t58+t60+t30*((t27*(t5-t9+t10-t36))/5.0-2.0)*pow(t5-t9+t10-t36,3.0)*(9.0/2.5E+1)+(t29*((t27*t42)/5.0+t28*(t5-t9+t10-t36)*(1.2E+1/5.0))*pow(t5-t9+t10-t36,3.0))/1.0E+2+t29*t38*t42*((t27*(t5-t9+t10-t36))/5.0-2.0)*(3.0/1.0E+2))-t2*(t31*t40*(3.0/1.25E+2)+t30*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*pow(t5-t9+t10-t36,3.0)*(1.8E+1/2.5E+1)+(t30*t42*pow(t5-t9+t10-t36,3.0))/5.0E+2+(t29*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*pow(t5-t9+t10-t36,3.0))/5.0E+1+t29*t38*t42*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*(3.0/5.0E+1))-(q0+Tqd0*(2.0/5.0)+TTqdd0/2.0E+1)*(t52+t53+t56+t58+t59+t60)+(t29*((t27*(t5-t9+t10-t36))/5.0-2.0)*pow(t5-t9+t10-t36,3.0))/1.0E+2-(t29*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*pow(t5-t9+t10-t36,3.0))/5.0E+1+t2*t31*t40*(3.0/1.25E+2)+t28*t38*t48*(3.0/1.0E+1)-q0*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*2.0E+1+(t2*t30*t42*pow(t5-t9+t10-t36,3.0))/5.0E+2;
}

double q_des_k_indep(double q0, double Tqd0, double TTqdd0, double kc, double s, double duration) {
    return q0 + Tqd0*s - 6*Tqd0*pow(s,3) + 8*Tqd0*pow(s,4) - 3*Tqd0*pow(s,5) + (TTqdd0*pow(s,2))*0.5 - (3*TTqdd0*pow(s,3))*0.5 + (3*TTqdd0*pow(s,4))*0.5 - (TTqdd0*pow(s,5))*0.5 + 10*kc*pow(s,3) - 15*kc*pow(s,4) + 6*kc*pow(s,5);
}

double qd_des_k_indep(double q0, double Tqd0, double TTqdd0, double kc, double s, double duration) {
    return (pow(s - 1,2)*(2*Tqd0 + 4*Tqd0*s + 2*TTqdd0*s - 30*Tqd0*pow(s,2) - 5*TTqdd0*pow(s,2) + 60*kc*pow(s,2)))*0.5 / duration;
}

double qdd_des_k_indep(double q0, double Tqd0, double TTqdd0, double kc, double s, double duration) {
    return -(s - 1.0)*(TTqdd0 - (36*Tqd0 + 8*TTqdd0)*s + (60*Tqd0 + 10*TTqdd0)*pow(s, 2) + 60*kc*s - 120*kc*pow(s,2))  / (duration * duration);
}

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR