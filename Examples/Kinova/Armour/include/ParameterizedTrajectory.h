#ifndef PARAMETERIZED_TRAJECTORY_H
#define PARAMETERIZED_TRAJECTORY_H

#include "RobotInfo.h"

namespace RAPTOR {
namespace Armour {

constexpr size_t NUM_TIME_STEPS = 128; // Number of time intervals partitioning the trajectory

// These values are specifically corresponded with a Bezuer curve parameterization
constexpr float QDD_DES_K_DEP_MAXIMA = (0.5 - sqrtf(3.0f) / 6);
constexpr float QDD_DES_K_DEP_MINIMA = (0.5 + sqrtf(3.0f) / 6);

// 5th order Bezier curve
// The initial position/velocity/acceleration is equal to q0/qd0/qdd0
// The end position is equal to q0 + k
// The end velocity/acceleration is equal to 0
//
// NOTE:
// This is just a simplified implementation!!!
// t is automatically set to range from 0 to 1, so you don't have to scale.
// Everything has to be changed if the range of t is not [0,1]
//
// -1 <= k <= 1
// k_actual = k * k_range
//
// q_des   = t^3*(6*t^2 - 15*t + 10) * k_actual + 
//           q0 + qd0*t - 6*qd0*t^3 + 8*qd0*t^4 - 3*qd0*t^5 + (qdd0*t^2)/2 - (3*qdd0*t^3)/2 + (3*qdd0*t^4)/2 - (qdd0*t^5)/2
//
// qd_des  = 30*t^2*(t - 1)^2 * k_actual + 
//           ((t - 1)^2*(2*qd0 + 4*qd0*t + 2*qdd0*t - 30*qd0*t^2 - 5*qdd0*t^2))/2
//
// qdd_des = 60*t*(2*t^2 - 3*t + 1) * k_actual + 
//           -(t - 1)*(qdd0 - 36*qd0*t - 8*qdd0*t + 60*qd0*t^2 + 10*qdd0*t^2)
//
class BezierCurveInterval{
public:
    using VecX = Eigen::VectorXf;
    using MatX = Eigen::MatrixXf;

    ultimate_bound ultimate_bound_info;

    VecX k_center;
    VecX k_range;

    float duration = 1.0;

    size_t num_time_steps = 0;
    std::vector<float> s_intervals;

    VecX q0;
    VecX qd0;
    VecX qdd0;

    VecX Tqd0; // qd0 * T
    VecX TTqdd0; // qdd0 * T ^ 2

    float q_des_k_indep_extrema_1[NUM_FACTORS] = {0.0};
    float q_des_k_indep_extrema_2[NUM_FACTORS] = {0.0};
    float q_des_k_indep_extremum_1[NUM_FACTORS] = {0.0};
    float q_des_k_indep_extremum_2[NUM_FACTORS] = {0.0};

    float qd_des_k_indep_extrema_1[NUM_FACTORS] = {0.0};
    float qd_des_k_indep_extrema_2[NUM_FACTORS] = {0.0};
    float qd_des_k_indep_extremum_1[NUM_FACTORS] = {0.0};
    float qd_des_k_indep_extremum_2[NUM_FACTORS] = {0.0};

    float qdd_des_k_indep_extrema_1[NUM_FACTORS] = {0.0};
    float qdd_des_k_indep_extrema_2[NUM_FACTORS] = {0.0};
    float qdd_des_k_indep_extremum_1[NUM_FACTORS] = {0.0};
    float qdd_des_k_indep_extremum_2[NUM_FACTORS] = {0.0};

    float ds = 0;

    // rotation matrix (and its transpose) of each joint
    PZsparseArray cos_q_des;
    PZsparseArray sin_q_des;
    PZsparseArray R;
    PZsparseArray R_t;

    // joint velocity
    PZsparseArray qd_des;

    // auxiliary joint velocity
    PZsparseArray qda_des;

    // joint acceleration
    PZsparseArray qdda_des;

    BezierCurveInterval();

    BezierCurveInterval(
        const VecX& q0_inp, 
        const VecX& qd0_inp, 
        const VecX& qdd0_inp,
        const VecX& k_center_inp,
        const VecX& k_range_inp,
        const float duration_inp,
        const ultimate_bound& ultimate_bound_info_inp);

    ~BezierCurveInterval() {};

    // set up the trajectory parameters
    void setTrajectoryParameters(
        const VecX& q0_inp, 
        const VecX& qd0_inp, 
        const VecX& qdd0_inp,
        const VecX& k_center_inp,
        const VecX& k_range_inp,
        const float duration_inp);

    // convert to polynomial zonotope using 1st/2nd order Taylor expansion
    void makePolyZono(const int s_ind);

    // return the min and max of the joint position throughout the whole desired trajectory
    void returnJointPositionExtremum(float* extremum, const float* k) const;

    // return the gradient of min and max of the joint position throughout the whole desired trajectory
    void returnJointPositionExtremumGradient(float* extremumGradient, const float* k) const;

    // return the min and max of the joint velocity throughout the whole desired trajectory
    void returnJointVelocityExtremum(float* extremum, const float* k) const;

    // return the gradient of min and max of the joint velocity throughout the whole desired trajectory
    void returnJointVelocityExtremumGradient(float* extremumGradient, const float* k) const;
};

// helper functions
// q0, qd0, qdd0, k here are scalars since all joints are using the same Bezier curve representation
float q_des_func(const float q0, const float Tqd0, const float TTqdd0, const float k, const float s);

float qd_des_func(const float q0, const float Tqd0, const float TTqdd0, const float k, const float s, const float T);

float qdd_des_func(const float q0, const float Tqd0, const float TTqdd0, const float k, const float s, const float T);

// derivative of the second extrema of q_des (when qd_des = 0) w.r.t k 
float q_des_extrema2_k_derivative(float q0, float Tqd0, float TTqdd0, float k);

// derivative of the third extrema of q_des (when qd_des = 0) w.r.t k 
float q_des_extrema3_k_derivative(float q0, float Tqd0, float TTqdd0, float k);

// derivative of the second extrema of qd_des (when qdd_des = 0) w.r.t k 
float qd_des_extrema2_k_derivative(float q0, float Tqd0, float TTqdd0, float k);

// derivative of the third extrema of qd_des (when qdd_des = 0) w.r.t k 
float qd_des_extrema3_k_derivative(float q0, float Tqd0, float TTqdd0, float k);

// k-independent part of q_des
float q_des_k_indep(float q0, float Tqd0, float TTqdd0, float kc, float s, float duration);

// k-independent part of qd_des
float qd_des_k_indep(float q0, float Tqd0, float TTqdd0, float kc, float s, float duration);

// k-independent part of qdd_des
float qdd_des_k_indep(float q0, float Tqd0, float TTqdd0, float kc, float s, float duration);

}; // namespace Armour
}; // namespace RAPTOR

#endif // PARAMETERIZED_TRAJECTORY_H