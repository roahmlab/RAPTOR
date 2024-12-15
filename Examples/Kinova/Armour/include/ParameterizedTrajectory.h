#ifndef PARAMETERIZED_TRAJECTORY_H
#define PARAMETERIZED_TRAJECTORY_H

#include "pinocchio/algorithm/model.hpp"
#include "pinocchio/algorithm/crba.hpp"

#include "RobotInfo.h"
#include "Utils.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

constexpr size_t DEFAULT_NUM_TIME_STEPS = 64; // Number of time intervals partitioning the trajectory

// These values are specifically corresponded with a Bezuer curve parameterization
constexpr double QDD_DES_K_DEP_MAXIMA = (0.5 - sqrtf(3.0) / 6);
constexpr double QDD_DES_K_DEP_MINIMA = (0.5 + sqrtf(3.0) / 6);

// 5th order Bezier curve
// The initial position/velocity/acceleration is equal to q0/q_d0/q_dd0
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
//           q0 + q_d0*t - 6*q_d0*t^3 + 8*q_d0*t^4 - 3*q_d0*t^5 + (q_dd0*t^2)/2 - (3*q_dd0*t^3)/2 + (3*q_dd0*t^4)/2 - (q_dd0*t^5)/2
//
// qd_des  = 30*t^2*(t - 1)^2 * k_actual + 
//           ((t - 1)^2*(2*q_d0 + 4*q_d0*t + 2*q_dd0*t - 30*q_d0*t^2 - 5*q_dd0*t^2))/2
//
// qdd_des = 60*t*(2*t^2 - 3*t + 1) * k_actual + 
//           -(t - 1)*(q_dd0 - 36*q_d0*t - 8*q_dd0*t + 60*q_d0*t^2 + 10*q_dd0*t^2)
//
class BezierCurveInterval{
public:
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    std::shared_ptr<RobotInfo> robotInfoPtr_ = nullptr;

    VecX k_center;
    VecX k_range;

    double duration = 1.0;

    size_t num_time_steps = 0;
    std::vector<double> s_intervals;

    VecX q0;
    VecX q_d0;
    VecX q_dd0;

    VecX Tqd0; // q_d0 * T
    VecX TTqdd0; // q_dd0 * T ^ 2

    double q_des_k_indep_extrema_1[NUM_FACTORS] = {0.0};
    double q_des_k_indep_extrema_2[NUM_FACTORS] = {0.0};
    double q_des_k_indep_extremum_1[NUM_FACTORS] = {0.0};
    double q_des_k_indep_extremum_2[NUM_FACTORS] = {0.0};

    double qd_des_k_indep_extrema_1[NUM_FACTORS] = {0.0};
    double qd_des_k_indep_extrema_2[NUM_FACTORS] = {0.0};
    double qd_des_k_indep_extremum_1[NUM_FACTORS] = {0.0};
    double qd_des_k_indep_extremum_2[NUM_FACTORS] = {0.0};

    double qdd_des_k_indep_extrema_1[NUM_FACTORS] = {0.0};
    double qdd_des_k_indep_extrema_2[NUM_FACTORS] = {0.0};
    double qdd_des_k_indep_extremum_1[NUM_FACTORS] = {0.0};
    double qdd_des_k_indep_extremum_2[NUM_FACTORS] = {0.0};

    double ds = 0;

    // rotation matrix (and its transpose) of each joint
    PZSparseMatrix q_des;

    // joint velocity
    PZSparseMatrix qd_des;

    // auxiliary joint velocity
    PZSparseMatrix qda_des;

    // joint acceleration
    PZSparseMatrix qdda_des;

    BezierCurveInterval();

    BezierCurveInterval(
        const VecX& q0_inp, 
        const VecX& q_d0_inp, 
        const VecX& q_dd0_inp,
        const VecX& k_center_inp,
        const VecX& k_range_inp,
        const double duration_inp,
        const std::shared_ptr<RobotInfo>& robotInfoPtr_inp,
        const size_t num_time_steps_inp = DEFAULT_NUM_TIME_STEPS);

    ~BezierCurveInterval() {};

    void sample_eigenvalues(size_t num_samples = 1000);

    // set up the trajectory parameters
    void setTrajectoryParameters(
        const VecX& q0_inp, 
        const VecX& q_d0_inp, 
        const VecX& q_dd0_inp,
        const VecX& k_center_inp,
        const VecX& k_range_inp,
        const double duration_inp);

    void computeTrajectories(
        const VecX& k,
        const double t,
        VecX& q,
        VecX& qd,
        VecX& qdd) const;

    // convert to polynomial zonotope using 1st/2nd order Taylor expansion
    void makePolyZono(const int s_ind);

    // return the min and max of the joint position throughout the whole desired trajectory
    void returnJointPositionExtremum(double* extremum, const double* k) const;

    // return the gradient of min and max of the joint position throughout the whole desired trajectory
    void returnJointPositionExtremumGradient(double* extremumGradient, const double* k) const;

    // return the min and max of the joint velocity throughout the whole desired trajectory
    void returnJointVelocityExtremum(double* extremum, const double* k) const;

    // return the gradient of min and max of the joint velocity throughout the whole desired trajectory
    void returnJointVelocityExtremumGradient(double* extremumGradient, const double* k) const;
};

// helper functions
// q0, q_d0, q_dd0, k here are scalars since all joints are using the same Bezier curve representation
double q_des_func(const double q0, const double Tqd0, const double TTqdd0, const double k, const double s);

double qd_des_func(const double q0, const double Tqd0, const double TTqdd0, const double k, const double s, const double T);

double qdd_des_func(const double q0, const double Tqd0, const double TTqdd0, const double k, const double s, const double T);

// derivative of the second extrema of q_des (when qd_des = 0) w.r.t k 
double q_des_extrema2_k_derivative(double q0, double Tqd0, double TTqdd0, double k);

// derivative of the third extrema of q_des (when qd_des = 0) w.r.t k 
double q_des_extrema3_k_derivative(double q0, double Tqd0, double TTqdd0, double k);

// derivative of the second extrema of qd_des (when qdd_des = 0) w.r.t k 
double qd_des_extrema2_k_derivative(double q0, double Tqd0, double TTqdd0, double k);

// derivative of the third extrema of qd_des (when qdd_des = 0) w.r.t k 
double qd_des_extrema3_k_derivative(double q0, double Tqd0, double TTqdd0, double k);

// k-independent part of q_des
double q_des_k_indep(double q0, double Tqd0, double TTqdd0, double kc, double s, double duration);

// k-independent part of qd_des
double qd_des_k_indep(double q0, double Tqd0, double TTqdd0, double kc, double s, double duration);

// k-independent part of qdd_des
double qdd_des_k_indep(double q0, double Tqd0, double TTqdd0, double kc, double s, double duration);

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR

#endif // PARAMETERIZED_TRAJECTORY_H