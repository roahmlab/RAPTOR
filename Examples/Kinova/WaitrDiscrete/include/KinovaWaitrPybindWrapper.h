#ifndef KINOVA_WAITR_PYBIND_WRAPPER_H
#define KINOVA_WAITR_PYBIND_WRAPPER_H

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include "KinovaWaitrOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

namespace RAPTOR {
namespace Kinova {

namespace nb = nanobind;

class KinovaWaitrPybindWrapper {
public:
    using Model = pinocchio::ModelTpl<float>;
    using Vec3 = Eigen::Vector3f;
    using VecX = Eigen::VectorXf;
    using MatX = Eigen::MatrixXf;

    using nb_1d_float = nb::ndarray<float, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
    using nb_2d_float = nb::ndarray<float, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

    // Constructor
    KinovaWaitrPybindWrapper() = default;

    KinovaWaitrPybindWrapper(const std::string urdf_filename,
                             const bool display_info);

    // Destructor
    ~KinovaWaitrPybindWrapper() = default;

    // Class methods
    void set_obstacles(const nb_2d_float obstacles_inp);

    void set_contact_surface_parameters(const float mu_inp,
                                        const float R_inp,
                                        const float maxSuctionForce_inp = 0.0,
                                        const float contactForceBuffer_inp = 0.0,
                                        const float frictionForceBuffer_inp = 0.0,
                                        const float ZMPBuffer_inp = 0.0);

    void set_end_effector(const nb_1d_float contact_position,
                          const float object_mass,
                          const nb_1d_float object_com,
                          const nb_2d_float object_inertia);

    void set_ipopt_parameters(const float tol,
                              const float constr_viol_tol,
                              const float obj_scaling_factor,
                              const float max_wall_time, 
                              const int print_level,
                              const std::string mu_strategy,
                              const std::string linear_solver,
                              const bool gradient_check);

    void set_trajectory_parameters(const nb_1d_float q0_inp,
                                   const nb_1d_float qd0_inp,
                                   const nb_1d_float qdd0_inp,
                                   const float duration_inp);

    void set_buffer(const nb_1d_float joint_limits_buffer_inp,
                    const nb_1d_float velocity_limits_buffer_inp,
                    const nb_1d_float torque_limits_buffer_inp);

    void set_target(const nb_1d_float q_des_inp,
                    const float tplan_inp);

    nb::tuple optimize();

    nb::tuple analyze_solution();

    // Class members
    // robot model
    Model model;

    int actual_model_nq = 0;

    // obstacle information
    int num_obstacles = 0;
    std::vector<Vec3> boxCenters;
    std::vector<Vec3> boxOrientation;
    std::vector<Vec3> boxSize;

    // trajectory information
    ArmourTrajectoryParameters atp;
    float T = 1;
    int N = 32;
    int degree = ARMOUR_BEZIER_CURVE_DEGREE;
    VecX qdes;
    float tplan = 0;
    int tplan_n = 0;

    // contact surface parameters
    circleContactSurfaceParams csp;

    // buffer information
    VecX joint_limits_buffer;
    VecX velocity_limits_buffer;
    VecX torque_limits_buffer;
    
    SmartPtr<KinovaWaitrOptimizer> mynlp;
    SmartPtr<IpoptApplication> app;

    // Flags to check if the parameters are set
    bool set_obstacles_check = false;
    bool set_contact_surface_parameters_check = false;
    bool set_end_effector_check = false;
    bool set_ipopt_parameters_check = false;
    bool set_trajectory_parameters_check = false;
    bool set_buffer_check = false;
    bool set_target_check = false;
    bool has_optimized = false;
};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // KINOVA_WAITR_PYBIND_WRAPPER_H
