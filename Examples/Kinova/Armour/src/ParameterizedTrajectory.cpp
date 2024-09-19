#include "ParameterizedTrajectory.h"

namespace RAPTOR {
namespace Armour {

BezierCurve::BezierCurve() {
    q0 = VecX::Zero(NUM_FACTORS);
    qd0 = VecX::Zero(NUM_FACTORS);
    qdd0 = VecX::Zero(NUM_FACTORS);
    Tqd0 = VecX::Zero(NUM_FACTORS);
    TTqdd0 = VecX::Zero(NUM_FACTORS);

    // compute a uniform time partition
    num_time_steps = NUM_TIME_STEPS;
    s_intervals.reserve(num_time_steps + 1);
    for (int i = 0; i <= num_time_steps; i++) {
        s_intervals.push_back(i / (float)num_time_steps);
    }

    // pre-allocate memory
    cos_q_des = PZsparseArray(NUM_FACTORS, num_time_steps);
    sin_q_des = PZsparseArray(NUM_FACTORS, num_time_steps);
    qd_des = PZsparseArray(NUM_FACTORS, num_time_steps);
    qda_des = PZsparseArray(NUM_FACTORS, num_time_steps);
    qdda_des = PZsparseArray(NUM_FACTORS, num_time_steps);

    ds = 1.0 / num_time_steps;
}

BezierCurve::BezierCurve(const VecX& q0_inp, 
                         const VecX& qd0_inp, 
                         const VecX& qdd0_inp,
                         const VecX& k_center_inp,
                         const VecX& k_range_inp,
                         const float duration_inp,
                         const ultimate_bound& ultimate_bound_info_inp) : 
    q0(q0_inp),
    qd0(qd0_inp),
    qdd0(qdd0_inp),
    k_center(k_center_inp),
    k_range(k_range_inp),
    duration(duration_inp),
    ultimate_bound_info(ultimate_bound_info_inp) {
    Tqd0 = qd0 * duration; 
    TTqdd0 = qdd0 * duration * duration; 

    if (false) {
        // compute a near-optimal time partition
        // first compute a finer trajectory for all joints at the center and record the extrema
        const float k[NUM_FACTORS] = {0.0};
        MatX qdd_traj = MatX::Zero(4, 2048);
        VecX qdd_min(4); qdd_min.setConstant(1e19);
        VecX qdd_max(4); qdd_max.setConstant(-1e19);
        for (int i = 0; i < 2048; i++) {
            float s = i / 2048.0;
            for (int j = 0; j < 4; j++) {
                qdd_traj(j, i) = qdd_des_func(q0[j], Tqd0[j], TTqdd0[j], k[j], s, duration);
            }
            qdd_min = qdd_min.cwiseMin(qdd_traj.col(i));
            qdd_max = qdd_max.cwiseMax(qdd_traj.col(i));
        }

        // compute the desired interval
        VecX qdd_range = qdd_max - qdd_min;
        float interval_size = qdd_range.maxCoeff() / 60;

        s_intervals.reserve(70 + 1);
        s_intervals.push_back(0.0);

        // forward propagate to compute the time partition
        int current_lb_idx = 0;
        int current_ub_idx = 0;
        float current_ub = 0.0;
        VecX current_qdd_min(NUM_FACTORS); current_qdd_min.setConstant(1e19);
        VecX current_qdd_max(NUM_FACTORS); current_qdd_max.setConstant(-1e19);
        while (current_ub_idx < 2048) {
            current_qdd_min = qdd_traj.col(current_lb_idx);
            current_qdd_max = qdd_traj.col(current_ub_idx);
            int p = current_lb_idx + 1;
            for (; p < 2048; p++) {
                current_qdd_min = current_qdd_min.cwiseMin(qdd_traj.col(p));
                current_qdd_max = current_qdd_max.cwiseMax(qdd_traj.col(p));
                if ((current_qdd_max - current_qdd_min).maxCoeff() > interval_size) {
                    current_ub_idx = p;
                    current_ub = current_ub_idx / 2048.0;
                    break;
                }
            }

            if (p == 2048) {
                s_intervals.push_back(1.0);
                break;
            }

            s_intervals.push_back(current_ub);

            current_lb_idx = current_ub_idx;
        }

        num_time_steps = s_intervals.size() - 1;
    }
    else {
        // compute a uniform time partition
        num_time_steps = NUM_TIME_STEPS;
        s_intervals.reserve(num_time_steps + 1);
        for (int i = 0; i <= num_time_steps; i++) {
            s_intervals.push_back(i / (float)num_time_steps);
        }
    }

    // pre-allocate memory
    cos_q_des = PZsparseArray(NUM_FACTORS, num_time_steps);
    sin_q_des = PZsparseArray(NUM_FACTORS, num_time_steps);
    qd_des = PZsparseArray(NUM_FACTORS, num_time_steps);
    qda_des = PZsparseArray(NUM_FACTORS, num_time_steps);
    qdda_des = PZsparseArray(NUM_FACTORS, num_time_steps);

    // initialize the extrema of the k independent part of q_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        q_des_k_indep_extrema_1[i] = (2*Tqd0[i] + TTqdd0[i] + sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 120*k_center[i]*Tqd0[i]))/(5*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        q_des_k_indep_extrema_2[i] = (2*Tqd0[i] + TTqdd0[i] - sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 120*k_center[i]*Tqd0[i]))/(5*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        q_des_k_indep_extremum_1[i] = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], q_des_k_indep_extrema_1[i], duration);
        q_des_k_indep_extremum_2[i] = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], q_des_k_indep_extrema_2[i], duration);
    }

    // initialize the extrema of the k independent part of qd_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        qd_des_k_indep_extrema_1[i] = (18*Tqd0[i] + 4*TTqdd0[i] + sqrt(6*(54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 20*TTqdd0[i]*k_center[i] - 180*Tqd0[i]*k_center[i] + 150*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qd_des_k_indep_extrema_2[i] = (18*Tqd0[i] + 4*TTqdd0[i] - sqrt(6*(54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 20*TTqdd0[i]*k_center[i] - 180*Tqd0[i]*k_center[i] + 150*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qd_des_k_indep_extremum_1[i] = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qd_des_k_indep_extrema_1[i], duration);
        qd_des_k_indep_extremum_2[i] = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qd_des_k_indep_extrema_2[i], duration);
    }

    // initialize the extrema of the k independent part of qdd_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        qdd_des_k_indep_extrema_1[i] = (32*Tqd0[i] + 6*TTqdd0[i] + sqrt(2*(152*pow(Tqd0[i],2) + 42*Tqd0[i]*TTqdd0[i] + 3*pow(TTqdd0[i],2) - 80*TTqdd0[i]*k_center[i] - 600*Tqd0[i]*k_center[i] + 600*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qdd_des_k_indep_extrema_2[i] = (32*Tqd0[i] + 6*TTqdd0[i] - sqrt(2*(152*pow(Tqd0[i],2) + 42*Tqd0[i]*TTqdd0[i] + 3*pow(TTqdd0[i],2) - 80*TTqdd0[i]*k_center[i] - 600*Tqd0[i]*k_center[i] + 600*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qdd_des_k_indep_extremum_1[i] = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qdd_des_k_indep_extrema_1[i], duration);
        qdd_des_k_indep_extremum_2[i] = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qdd_des_k_indep_extrema_2[i], duration);
    }

    ds = 1.0 / num_time_steps;
}

void BezierCurve::setTrajectoryParameters(const VecX& q0_inp, 
                                          const VecX& qd0_inp, 
                                          const VecX& qdd0_inp,
                                          const VecX& k_center_inp,
                                          const VecX& k_range_inp,
                                          const float duration_inp) {
    q0 = q0_inp;
    qd0 = qd0_inp;
    qdd0 = qdd0_inp;
    k_center = k_center_inp;
    k_range = k_range_inp;
    duration = duration_inp;

    Tqd0 = qd0 * duration; 
    TTqdd0 = qdd0 * duration * duration;   

    // initialize the extrema of the k independent part of q_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        q_des_k_indep_extrema_1[i] = (2*Tqd0[i] + TTqdd0[i] + sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 120*k_center[i]*Tqd0[i]))/(5*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        q_des_k_indep_extrema_2[i] = (2*Tqd0[i] + TTqdd0[i] - sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 120*k_center[i]*Tqd0[i]))/(5*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        q_des_k_indep_extremum_1[i] = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], q_des_k_indep_extrema_1[i], duration);
        q_des_k_indep_extremum_2[i] = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], q_des_k_indep_extrema_2[i], duration);
    }

    // initialize the extrema of the k independent part of qd_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        qd_des_k_indep_extrema_1[i] = (18*Tqd0[i] + 4*TTqdd0[i] + sqrt(6*(54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 20*TTqdd0[i]*k_center[i] - 180*Tqd0[i]*k_center[i] + 150*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qd_des_k_indep_extrema_2[i] = (18*Tqd0[i] + 4*TTqdd0[i] - sqrt(6*(54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2) - 20*TTqdd0[i]*k_center[i] - 180*Tqd0[i]*k_center[i] + 150*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qd_des_k_indep_extremum_1[i] = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qd_des_k_indep_extrema_1[i], duration);
        qd_des_k_indep_extremum_2[i] = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qd_des_k_indep_extrema_2[i], duration);
    }

    // initialize the extrema of the k independent part of qdd_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        qdd_des_k_indep_extrema_1[i] = (32*Tqd0[i] + 6*TTqdd0[i] + sqrt(2*(152*pow(Tqd0[i],2) + 42*Tqd0[i]*TTqdd0[i] + 3*pow(TTqdd0[i],2) - 80*TTqdd0[i]*k_center[i] - 600*Tqd0[i]*k_center[i] + 600*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qdd_des_k_indep_extrema_2[i] = (32*Tqd0[i] + 6*TTqdd0[i] - sqrt(2*(152*pow(Tqd0[i],2) + 42*Tqd0[i]*TTqdd0[i] + 3*pow(TTqdd0[i],2) - 80*TTqdd0[i]*k_center[i] - 600*Tqd0[i]*k_center[i] + 600*pow(k_center[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i] - 12*k_center[i]));
        qdd_des_k_indep_extremum_1[i] = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qdd_des_k_indep_extrema_1[i], duration);
        qdd_des_k_indep_extremum_2[i] = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], qdd_des_k_indep_extrema_2[i], duration);
    }                                          
}

void BezierCurve::makePolyZono(const int s_ind) {
    assert(s_ind < num_time_steps);

    // const float s_lb = s_ind * ds;
    // const float s_ub = (s_ind + 1) * ds;
    const float s_lb = s_intervals[s_ind];
    const float s_ub = s_intervals[s_ind + 1];

    const Interval t_int(s_lb, s_ub);

    for (int i = 0; i < NUM_FACTORS; i++) {
        const float k_range_elt = k_range(i);

        // Part 1: q_des
        float k_dep_coeff_lb = pow(s_lb,3) * (6 * pow(s_lb,2) - 15 * s_lb + 10);
        float k_dep_coeff_ub = pow(s_ub,3) * (6 * pow(s_ub,2) - 15 * s_ub + 10);
        float k_dep_coeff_center = (k_dep_coeff_ub + k_dep_coeff_lb) * 0.5;
        float k_dep_coeff_radius = (k_dep_coeff_ub - k_dep_coeff_lb) * 0.5 * k_range_elt;
        
        float k_indep_lb = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], s_lb, duration);
        float k_indep_ub = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], k_center[i], s_ub, duration);
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
        float k_indep_radius = (k_indep_ub - k_indep_lb) * 0.5;
        float q_des_center = (k_indep_lb + k_indep_ub) * 0.5;
        
        // q_des_k_dep = k_dep_coeff_center * k;
        Interval q_des_radius_int(-k_dep_coeff_radius - k_indep_radius - ultimate_bound_info.qe, 
                                  k_dep_coeff_radius + k_indep_radius + ultimate_bound_info.qe);
        
        // q_des_int = q_des_center + q_des_k_dep + q_des_radius_int;

        // first order Taylor expansion
        // Part 1.a: cosf(q_des) 
        float cos_q_des_center = cosf(q_des_center);
        Interval cos_q_des_radius_int = - q_des_radius_int * sinf(q_des_center) 
                                        - 0.5f * cos(q_des_center + k_dep_coeff_center * Interval(-k_range_elt, k_range_elt) + q_des_radius_int) 
                                            * pow(q_des_radius_int + k_dep_coeff_center * Interval(-k_range_elt, k_range_elt), 2);

        cos_q_des_center += getCenter(cos_q_des_radius_int);
        cos_q_des_radius_int = cos_q_des_radius_int - getCenter(cos_q_des_radius_int);
        float cos_q_des_coeff[] = {-k_dep_coeff_center * k_range_elt * sinf(q_des_center), getRadius(cos_q_des_radius_int)}; 

        // cos_q_des_int = cos_q_des_center + cos_q_des_coeff[0] * k + cos_q_des_coeff[1] * cosqe;
        uint32_t cos_q_des_degree[2][NUM_FACTORS * 6] = {0};
        cos_q_des_degree[0][i] = 1; // k
        cos_q_des_degree[1][i + NUM_FACTORS * 4] = 1; // cosqe

        cos_q_des(i, s_ind) = PZsparse(cos_q_des_center, cos_q_des_coeff, cos_q_des_degree, 2);

        // Part 1.b: sinf(q_des) 
        float sin_q_des_center = sinf(q_des_center);
        Interval sin_q_des_radius_int = q_des_radius_int * cosf(q_des_center) 
                                        - 0.5f * sin(q_des_center + k_dep_coeff_center * Interval(-k_range_elt, k_range_elt) + q_des_radius_int) 
                                            * pow(q_des_radius_int + k_dep_coeff_center * Interval(-k_range_elt, k_range_elt), 2);

        sin_q_des_center += getCenter(sin_q_des_radius_int);
        sin_q_des_radius_int = sin_q_des_radius_int - getCenter(sin_q_des_radius_int);
        float sin_q_des_coeff[] = {k_dep_coeff_center * k_range_elt * cosf(q_des_center), getRadius(sin_q_des_radius_int)};

        // sin_q_des_int = sin_q_des_center + sin_q_des_coeff[0] * k + sin_q_des_coeff[1] * sinqe;
        uint32_t sin_q_des_degree[2][NUM_FACTORS * 6] = {0};
        sin_q_des_degree[0][i] = 1; // k
        sin_q_des_degree[1][i + NUM_FACTORS * 5] = 1; // sinqe

        sin_q_des(i, s_ind) = PZsparse(sin_q_des_center, sin_q_des_coeff, sin_q_des_degree, 2);

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

        k_dep_coeff_center = (k_dep_coeff_ub + k_dep_coeff_lb) * 0.5 * k_range_elt; // Have to scale to [-1,1] in order to fit in PZsparse
        k_dep_coeff_radius = (k_dep_coeff_ub - k_dep_coeff_lb) * 0.5 * k_range_elt; // Have to scale to [-1,1] in order to fit in PZsparse

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
        float qd_des_center = (k_indep_lb + k_indep_ub) * 0.5;

        float qd_des_coeff[] = {k_dep_coeff_center, k_dep_coeff_radius + k_indep_radius + ultimate_bound_info.qde};

        uint32_t qd_des_degree[2][NUM_FACTORS * 6] = {0};
        qd_des_degree[0][i] = 1; // k
        qd_des_degree[1][i + NUM_FACTORS * 1] = 1; // qde

        // qd_des_int = qd_des_center + qd_des_coeff[0] * k + qd_des_coeff[1] * qde;
        qd_des(i, s_ind) = PZsparse(qd_des_center, qd_des_coeff, qd_des_degree, 2);

        float qda_des_coeff[] = {k_dep_coeff_center, k_dep_coeff_radius + k_indep_radius + ultimate_bound_info.qdae};

        uint32_t qda_des_degree[2][NUM_FACTORS * 6] = {0};
        qda_des_degree[0][i] = 1; // k
        qda_des_degree[1][i + NUM_FACTORS * 2] = 1; // qdae

        // qda_des_int = qd_des_center + qda_des_coeff[0] * k + qda_des_coeff[1] * qdae;
        qda_des(i, s_ind) = PZsparse(qd_des_center, qda_des_coeff, qda_des_degree, 2);

        // Part 3: qdd_des
        float temp_lb = (60 * s_lb * (2 * pow(s_lb,2) - 3 * s_lb + 1)) / duration / duration;
        float temp_ub = (60 * s_ub * (2 * pow(s_ub,2) - 3 * s_ub + 1)) / duration / duration;
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
        
        k_dep_coeff_center = (k_dep_coeff_ub + k_dep_coeff_lb) * 0.5 * k_range_elt; // Have to scale to [-1,1] in order to fit in PZsparse
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
        float qdd_des_center = (k_indep_lb + k_indep_ub) * 0.5;

        float qdd_des_coeff[] = {k_dep_coeff_center, k_dep_coeff_radius + k_indep_radius + ultimate_bound_info.qddae};

        uint32_t qdd_des_degree[2][NUM_FACTORS * 6] = {0};
        qdd_des_degree[0][i] = 1; // k
        qdd_des_degree[1][i + NUM_FACTORS * 3] = 1; // qddae

        // qdd_des_int = qdd_des_center + qdd_des_coeff[0] * k + qdd_des_coeff[1] * qdde;
        qdda_des(i, s_ind) = PZsparse(qdd_des_center, qdd_des_coeff, qdd_des_degree, 2);
    }
}

void BezierCurve::returnJointPositionExtremum(float* extremum, const float* k) const {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        float k_actual = k_center(i) + k_range(i) * k[i];

        // list all possible extremas
        float extrema1 = 0;
        float extrema2 = (2*Tqd0[i] + TTqdd0[i] + sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] - 120*k_actual*Tqd0[i] + pow(TTqdd0[i],2))) / (5*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        float extrema3 = (2*Tqd0[i] + TTqdd0[i] - sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] - 120*k_actual*Tqd0[i] + pow(TTqdd0[i],2))) / (5*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        float extrema4 = 1;

        // get extremums of all extremas
        float extremum1 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema1);
        float extremum2 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema2);
        float extremum3 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema3);
        float extremum4 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema4);

        // find the min and max values
        float minPosition = std::min(extremum1, extremum4);
        float maxPosition = std::max(extremum1, extremum4);
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

void BezierCurve::returnJointPositionExtremumGradient(float* extremumGradient, const float* k) const {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        float k_actual = k_center(i) + k_range(i) * k[i];

        // list all possible extremas
        float extrema1 = 0;
        float extrema2 = (2*Tqd0[i] + TTqdd0[i] + sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] - 120*k_actual*Tqd0[i] + pow(TTqdd0[i],2))) / (5*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        float extrema3 = (2*Tqd0[i] + TTqdd0[i] - sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] - 120*k_actual*Tqd0[i] + pow(TTqdd0[i],2))) / (5*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        float extrema4 = 1;

        // get extremums of all extremas
        float extremum1 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema1);
        float extremum2 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema2);
        float extremum3 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema3);
        float extremum4 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema4);

        // find the min and max values
        float minPosition;
        int minId;
        float maxPosition;
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

        float minPositionGradient;
        float maxPositionGradient;

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

void BezierCurve::returnJointVelocityExtremum(float* extremum, const float* k) const {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        float k_actual = k_center(i) + k_range(i) * k[i];

        // list all possible extremas
        float extrema1 = 0;
        float extrema2 = (18*Tqd0[i] - 30*k_actual + 4*TTqdd0[i] + sqrt(6*(150*pow(k_actual,2) - 180*k_actual*Tqd0[i] - 20*k_actual*TTqdd0[i] + 54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        float extrema3 = (18*Tqd0[i] - 30*k_actual + 4*TTqdd0[i] - sqrt(6*(150*pow(k_actual,2) - 180*k_actual*Tqd0[i] - 20*k_actual*TTqdd0[i] + 54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        float extrema4 = 1;

        // get extremums of all extremas
        float extremum1 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema1, duration);
        float extremum2 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema2, duration);
        float extremum3 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema3, duration);
        float extremum4 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema4, duration);

        // find the min and max values
        float minVelocity = std::min(extremum1, extremum4);
        float maxVelocity = std::max(extremum1, extremum4);
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

void BezierCurve::returnJointVelocityExtremumGradient(float* extremumGradient, const float* k) const {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        float k_actual = k_center(i) + k_range(i) * k[i];

        // list all possible extremas
        float extrema1 = 0;
        float extrema2 = (18*Tqd0[i] - 30*k_actual + 4*TTqdd0[i] + sqrt(6*(150*pow(k_actual,2) - 180*k_actual*Tqd0[i] - 20*k_actual*TTqdd0[i] + 54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        float extrema3 = (18*Tqd0[i] - 30*k_actual + 4*TTqdd0[i] - sqrt(6*(150*pow(k_actual,2) - 180*k_actual*Tqd0[i] - 20*k_actual*TTqdd0[i] + 54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        float extrema4 = 1;

        // get extremums of all extremas
        float extremum1 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema1, duration);
        float extremum2 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema2, duration);
        float extremum3 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema3, duration);
        float extremum4 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema4, duration);

        // find the min and max values
        float minVelocity;
        int minId;
        float maxVelocity;
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

        float minVelocityGradient;
        float maxVelocityGradient;

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
                minVelocityGradient = 1.0;
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
                maxVelocityGradient = 1.0;
                break;
            default:
                break;
        }

        for (int j = 0; j < NUM_FACTORS; j++) {
            if (i == j) {
                extremumGradient[(i              ) * NUM_FACTORS + j] = minVelocityGradient * k_range(i);
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = maxVelocityGradient * k_range(i);
            }
            else {
                extremumGradient[(i              ) * NUM_FACTORS + j] = 0.0;
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = 0.0;
            }
        }
    }
}

float q_des_func(const float q0, const float Tqd0, const float TTqdd0, const float k, const float s) {
    float B0 = -pow(s - 1,5);
    float B1 = 5*s*pow(s - 1,4);
    float B2 = -10*pow(s,2)*pow(s - 1,3);
    float B3 = 10*pow(s,3)*pow(s - 1,2);
    float B4 = -5*pow(s,4)*(s - 1);
    float B5 = pow(s,5);
    float beta0 = q0;
    float beta1 = q0 + Tqd0/5;
    float beta2 = q0 + (2*Tqd0)/5 + TTqdd0/20;
    float beta3 = q0 + k;
    float beta4 = q0 + k;
    float beta5 = q0 + k;
    return B0 * beta0 + B1 * beta1 + B2 * beta2 + B3 * beta3 + B4 * beta4 + B5 * beta5;
}

float qd_des_func(const float q0, const float Tqd0, const float TTqdd0, const float k, const float s, const float T) {
    float dB0 = pow(s-1.0,4.0)*-5.0;
    float dB1 = s*pow(s-1.0,3.0)*2.0E+1+pow(s-1.0,4.0)*5.0;
    float dB2 = s*pow(s-1.0,3.0)*-2.0E+1-(s*s)*pow(s-1.0,2.0)*3.0E+1;
    float dB3 = pow(s,3.0)*(s*2.0-2.0)*1.0E+1+(s*s)*pow(s-1.0,2.0)*3.0E+1;
    float dB4 = pow(s,3.0)*(s-1.0)*-2.0E+1-pow(s,4.0)*5.0;
    float dB5 = pow(s,4.0)*5.0;
    float beta0 = q0;
    float beta1 = q0 + Tqd0/5;
    float beta2 = q0 + (2*Tqd0)/5 + TTqdd0/20;
    float beta3 = q0 + k;
    float beta4 = q0 + k;
    float beta5 = q0 + k;
    float traj = dB0 * beta0 + dB1 * beta1 + dB2 * beta2 + dB3 * beta3 + dB4 * beta4 + dB5 * beta5;
    return traj / T;
}

float qdd_des_func(const float q0, const float Tqd0, const float TTqdd0, const float k, const float s, const float T) {
    float t2 = s*2.0;
    float t3 = s*s;
    float t4 = s*s*s;
    float t5 = s-1.0;
    float t6 = t2-2.0;
    float t7 = t4*2.0E+1;
    float t8 = t5*t5;
    float t9 = t5*t5*t5;
    float t10 = t9*2.0E+1;
    float t11 = s*t8*6.0E+1;
    float t12 = -t10;
    float ddB0 = t12;
    float ddB1 = t9*4.0E+1+t11;
    float ddB2 = t12-s*t8*1.2E+2-t3*t6*3.0E+1;
    float ddB3 = t7+t11+t3*t6*6.0E+1;
    float ddB4 = t4*-4.0E+1-t3*t5*6.0E+1;
    float ddB5 = t7;
    float beta0 = q0;
    float beta1 = q0 + Tqd0/5;
    float beta2 = q0 + (2*Tqd0)/5 + TTqdd0/20;
    float beta3 = q0 + k;
    float beta4 = q0 + k;
    float beta5 = q0 + k;
    float traj = ddB0 * beta0 + ddB1 * beta1 + ddB2 * beta2 + ddB3 * beta3 + ddB4 * beta4 + ddB5 * beta5;
    return traj / (T * T);
}

float q_des_extrema2_k_derivative(float q0, float Tqd0, float TTqdd0, float k) {
    float t2 = k+q0;
    float t3 = Tqd0*2.0;
    float t4 = Tqd0*6.0;
    float t5 = Tqd0*Tqd0;
    float t6 = TTqdd0*TTqdd0;
    float t7 = k*1.2E+1;
    float t8 = Tqd0*TTqdd0*1.4E+1;
    float t10 = Tqd0/5.0;
    float t11 = Tqd0*(2.0/5.0);
    float t12 = k*Tqd0*1.2E+2;
    float t13 = TTqdd0/2.0E+1;
    float t9 = -t7;
    float t14 = t5*6.4E+1;
    float t15 = -t12;
    float t16 = q0+t10;
    float t18 = q0+t11+t13;
    float t17 = TTqdd0+t4+t9;
    float t24 = t6+t8+t14+t15;
    float t19 = 1.0/t17;
    float t25 = sqrt(t24);
    float t20 = t19*t19;
    float t21 = t19*t19*t19;
    float t23 = t19*t19*t19*t19*t19;
    float t26 = 1.0/t25;
    float t27 = TTqdd0+t3+t25;
    float t22 = t20*t20;
    float t28 = t27*t27;
    float t29 = t27*t27*t27;
    float t31 = t27*t27*t27*t27*t27;
    float t32 = Tqd0*t19*t26*1.2E+1;
    float t34 = (t19*t27)/5.0;
    float t35 = t20*t27*(1.2E+1/5.0);
    float t30 = t28*t28;
    float t33 = -t32;
    float t36 = t34-1.0;
    float t37 = t36*t36;
    float t38 = t36*t36*t36;
    float t40 = t33+t35;
    float t39 = t37*t37;
    return (t23*t31)/3.125E+3+t2*(t20*t20*t20)*t31*(1.2E+1/6.25E+2)+t21*t29*t37*(2.0/2.5E+1)-(t22*t30*t36)/1.25E+2+q0*t39*(t32-t35)*5.0+t2*t22*t29*t37*(7.2E+1/2.5E+1)-t2*t23*t30*t36*(4.8E+1/1.25E+2)+t16*t20*t27*t39*1.2E+1-t18*t21*t28*t38*(4.8E+1/5.0)+(t2*t22*t30*(t32-t35))/1.25E+2-Tqd0*t2*t23*t26*t30*(1.2E+1/1.25E+2)-Tqd0*t16*t19*t26*t39*6.0E+1-t2*t21*t29*t36*(t32-t35)*(4.0/2.5E+1)-t16*t19*t27*t38*(t32-t35)*4.0+t18*t20*t28*t37*(t32-t35)*(6.0/5.0)-Tqd0*t2*t21*t26*t28*t37*(7.2E+1/5.0)+Tqd0*t2*t22*t26*t29*t36*(4.8E+1/2.5E+1)+Tqd0*t18*t20*t26*t27*t38*4.8E+1;
}

float q_des_extrema3_k_derivative(float q0, float Tqd0, float TTqdd0, float k) {
    float t2 = k+q0;
    float t3 = Tqd0*2.0;
    float t4 = Tqd0*6.0;
    float t5 = Tqd0*Tqd0;
    float t6 = TTqdd0*TTqdd0;
    float t7 = k*1.2E+1;
    float t8 = Tqd0*TTqdd0*1.4E+1;
    float t10 = Tqd0/5.0;
    float t11 = Tqd0*(2.0/5.0);
    float t12 = k*Tqd0*1.2E+2;
    float t13 = TTqdd0/2.0E+1;
    float t9 = -t7;
    float t14 = t5*6.4E+1;
    float t15 = -t12;
    float t16 = q0+t10;
    float t18 = q0+t11+t13;
    float t17 = TTqdd0+t4+t9;
    float t24 = t6+t8+t14+t15;
    float t19 = 1.0/t17;
    float t25 = sqrt(t24);
    float t20 = t19*t19;
    float t21 = t19*t19*t19;
    float t23 = t19*t19*t19*t19*t19;
    float t26 = 1.0/t25;
    float t27 = -t25;
    float t22 = t20*t20;
    float t28 = TTqdd0+t3+t27;
    float t33 = Tqd0*t19*t26*1.2E+1;
    float t29 = t28*t28;
    float t30 = t28*t28*t28;
    float t32 = t28*t28*t28*t28*t28;
    float t34 = (t19*t28)/5.0;
    float t35 = t20*t28*(1.2E+1/5.0);
    float t31 = t29*t29;
    float t36 = t34-1.0;
    float t40 = t33+t35;
    float t37 = t36*t36;
    float t38 = t36*t36*t36;
    float t39 = t37*t37;
    return (t23*t32)/3.125E+3+t2*(t20*t20*t20)*t32*(1.2E+1/6.25E+2)-q0*t39*t40*5.0+t21*t30*t37*(2.0/2.5E+1)-(t22*t31*t36)/1.25E+2+t2*t22*t30*t37*(7.2E+1/2.5E+1)-t2*t23*t31*t36*(4.8E+1/1.25E+2)-(t2*t22*t31*t40)/1.25E+2+t16*t20*t28*t39*1.2E+1-t18*t21*t29*t38*(4.8E+1/5.0)+Tqd0*t2*t23*t26*t31*(1.2E+1/1.25E+2)+Tqd0*t16*t19*t26*t39*6.0E+1+t2*t21*t30*t36*t40*(4.0/2.5E+1)+t16*t19*t28*t38*t40*4.0-t18*t20*t29*t37*t40*(6.0/5.0)+Tqd0*t2*t21*t26*t29*t37*(7.2E+1/5.0)-Tqd0*t2*t22*t26*t30*t36*(4.8E+1/2.5E+1)-Tqd0*t18*t20*t26*t28*t38*4.8E+1;
}

float qd_des_extrema2_k_derivative(float q0, float Tqd0, float TTqdd0, float k) {
    float t2 = k+q0;
    float t3 = k*k;
    float t4 = Tqd0*6.0;
    float t5 = TTqdd0*4.0;
    float t6 = Tqd0*Tqd0;
    float t7 = TTqdd0*TTqdd0;
    float t8 = k*1.2E+1;
    float t9 = k*3.0E+1;
    float t10 = Tqd0*1.8E+1;
    float t11 = TTqdd0*2.0E+1;
    float t13 = Tqd0*TTqdd0*1.4E+1;
    float t14 = sqrt(6.0);
    float t17 = k*3.0E+2;
    float t18 = Tqd0*1.8E+2;
    float t21 = k*TTqdd0*-2.0E+1;
    float t24 = k*Tqd0*-1.8E+2;
    float t12 = k*t11;
    float t15 = -t8;
    float t16 = -t9;
    float t19 = t6*5.4E+1;
    float t20 = k*t18;
    float t22 = -t17;
    float t23 = t3*1.5E+2;
    float t25 = TTqdd0+t4+t15;
    float t31 = t11+t18+t22;
    float t32 = t7+t13+t19+t21+t23+t24;
    float t26 = 1.0/t25;
    float t33 = sqrt(t32);
    float t27 = t26*t26;
    float t28 = t26*t26*t26;
    float t30 = t26*t26*t26*t26*t26;
    float t34 = 1.0/t33;
    float t35 = t14*t33;
    float t29 = t27*t27;
    float t36 = t5+t10+t16+t35;
    float t40 = (t14*t31*t34)/2.0;
    float t37 = t36*t36;
    float t38 = t36*t36*t36;
    float t41 = t40+3.0E+1;
    float t42 = (t26*t36)/5.0;
    float t43 = t27*t36*(6.0/5.0);
    float t44 = (t26*t36)/1.0E+1;
    float t39 = t37*t37;
    float t45 = t42-2.0;
    float t46 = t44-1.0;
    float t49 = (t26*t41)/1.0E+1;
    float t47 = t46*t46;
    float t48 = t46*t46*t46;
    float t50 = -t49;
    float t51 = t27*t36*t48*2.4E+1;
    float t52 = t28*t37*t47*(3.6E+1/5.0);
    float t53 = t43+t50;
    float t54 = t26*t41*t48*2.0;
    float t56 = t27*t36*t41*t47*(3.0/5.0);
    float t55 = -t54;
    float t57 = -t56;
    float t58 = t26*t36*t47*t53*6.0;
    float t59 = t27*t37*t46*t53*(3.0/5.0);
    return (q0+Tqd0/5.0)*(t51+t55+t58+t48*t53*2.0E+1)+t2*(t52+t57+t59+(t28*t38*(t27*t36*(1.2E+1/5.0)-(t26*t41)/5.0))/1.0E+2+t29*t38*t45*(9.0/2.5E+1)-t28*t37*t41*t45*(3.0/1.0E+2))-t2*(t30*t39*(3.0/1.25E+2)-(t29*t38*t41)/5.0E+2+t29*t38*t46*(1.8E+1/2.5E+1)+(t28*t38*t53)/5.0E+1-t28*t37*t41*t46*(3.0/5.0E+1))-(q0+Tqd0*(2.0/5.0)+TTqdd0/2.0E+1)*(t51+t52+t55+t57+t58+t59)-q0*t48*t53*2.0E+1+t2*t30*t39*(3.0/1.25E+2)+t27*t37*t47*(3.0/1.0E+1)+(t28*t38*t45)/1.0E+2-(t28*t38*t46)/5.0E+1-(t2*t29*t38*t41)/5.0E+2;
}

float qd_des_extrema3_k_derivative(float q0, float Tqd0, float TTqdd0, float k) {
    float t2 = k+q0;
    float t3 = k*k;
    float t4 = Tqd0*6.0;
    float t5 = TTqdd0*4.0;
    float t6 = Tqd0*Tqd0;
    float t7 = TTqdd0*TTqdd0;
    float t8 = k*1.2E+1;
    float t9 = k*3.0E+1;
    float t10 = Tqd0*1.8E+1;
    float t12 = TTqdd0*2.0E+1;
    float t14 = Tqd0*TTqdd0*1.4E+1;
    float t15 = sqrt(6.0);
    float t17 = k*3.0E+2;
    float t19 = Tqd0*1.8E+2;
    float t22 = k*TTqdd0*-2.0E+1;
    float t25 = k*Tqd0*-1.8E+2;
    float t11 = -t5;
    float t13 = k*t12;
    float t16 = -t8;
    float t18 = -t10;
    float t20 = t6*5.4E+1;
    float t21 = k*t19;
    float t23 = -t17;
    float t24 = t3*1.5E+2;
    float t26 = TTqdd0+t4+t16;
    float t32 = t12+t19+t23;
    float t33 = t7+t14+t20+t22+t24+t25;
    float t27 = 1.0/t26;
    float t34 = sqrt(t33);
    float t28 = t27*t27;
    float t29 = t27*t27*t27;
    float t31 = t27*t27*t27*t27*t27;
    float t35 = 1.0/t34;
    float t36 = t15*t34;
    float t30 = t28*t28;
    float t37 = t9+t11+t18+t36;
    float t38 = pow(t5-t9+t10-t36,2.0);
    float t39 = -pow(t5-t9+t10-t36,3.0);
    float t41 = (t15*t32*t35)/2.0;
    float t43 = t27*(t5-t9+t10-t36)*(-1.0/5.0);
    float t44 = t28*(t5-t9+t10-t36)*(-6.0/5.0);
    float t45 = t27*(t5-t9+t10-t36)*(-1.0/1.0E+1);
    float t48 = pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,2.0);
    float t49 = -pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0);
    float t52 = t28*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*(t5-t9+t10-t36)*2.4E+1;
    float t40 = t38*t38;
    float t42 = t41-3.0E+1;
    float t46 = t43+2.0;
    float t47 = t45+1.0;
    float t53 = t29*t38*t48*(3.6E+1/5.0);
    float t50 = (t27*t42)/1.0E+1;
    float t55 = t27*t42*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*-2.0;
    float t56 = t27*t42*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*2.0;
    float t57 = t28*t42*t48*(t5-t9+t10-t36)*(-3.0/5.0);
    float t58 = t28*t42*t48*(t5-t9+t10-t36)*(3.0/5.0);
    float t51 = -t50;
    float t59 = t27*t48*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*(t5-t9+t10-t36)*6.0;
    float t60 = t28*t38*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*(3.0/5.0);
    float t54 = t44+t51;
    return (q0+Tqd0/5.0)*(t52+t56+t59+pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*2.0E+1)+t2*(t53+t58+t60+t30*((t27*(t5-t9+t10-t36))/5.0-2.0)*pow(t5-t9+t10-t36,3.0)*(9.0/2.5E+1)+(t29*((t27*t42)/5.0+t28*(t5-t9+t10-t36)*(1.2E+1/5.0))*pow(t5-t9+t10-t36,3.0))/1.0E+2+t29*t38*t42*((t27*(t5-t9+t10-t36))/5.0-2.0)*(3.0/1.0E+2))-t2*(t31*t40*(3.0/1.25E+2)+t30*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*pow(t5-t9+t10-t36,3.0)*(1.8E+1/2.5E+1)+(t30*t42*pow(t5-t9+t10-t36,3.0))/5.0E+2+(t29*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*pow(t5-t9+t10-t36,3.0))/5.0E+1+t29*t38*t42*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*(3.0/5.0E+1))-(q0+Tqd0*(2.0/5.0)+TTqdd0/2.0E+1)*(t52+t53+t56+t58+t59+t60)+(t29*((t27*(t5-t9+t10-t36))/5.0-2.0)*pow(t5-t9+t10-t36,3.0))/1.0E+2-(t29*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*pow(t5-t9+t10-t36,3.0))/5.0E+1+t2*t31*t40*(3.0/1.25E+2)+t28*t38*t48*(3.0/1.0E+1)-q0*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*2.0E+1+(t2*t30*t42*pow(t5-t9+t10-t36,3.0))/5.0E+2;
}

float q_des_k_indep(float q0, float Tqd0, float TTqdd0, float kc, float s, float duration) {
    return q0 + Tqd0*s - 6*Tqd0*pow(s,3) + 8*Tqd0*pow(s,4) - 3*Tqd0*pow(s,5) + (TTqdd0*pow(s,2))*0.5 - (3*TTqdd0*pow(s,3))*0.5 + (3*TTqdd0*pow(s,4))*0.5 - (TTqdd0*pow(s,5))*0.5 + 10*kc*pow(s,3) - 15*kc*pow(s,4) + 6*kc*pow(s,5);
}

float qd_des_k_indep(float q0, float Tqd0, float TTqdd0, float kc, float s, float duration) {
    return (pow(s - 1,2)*(2*Tqd0 + 4*Tqd0*s + 2*TTqdd0*s - 30*Tqd0*pow(s,2) - 5*TTqdd0*pow(s,2) + 60*kc*pow(s,2)))*0.5 / duration;
}

float qdd_des_k_indep(float q0, float Tqd0, float TTqdd0, float kc, float s, float duration) {
    return -(s - 1.0)*(TTqdd0 - (36*Tqd0 + 8*TTqdd0)*s + (60*Tqd0 + 10*TTqdd0)*pow(s, 2) + 60*kc*s - 120*kc*pow(s,2))  / (duration * duration);
}

}; // namespace Armour
}; // namespace RAPTOR