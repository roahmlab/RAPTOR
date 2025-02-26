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
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    using nb_1d_double = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
    using nb_2d_double = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

    // Constructor
    KinovaWaitrPybindWrapper() = default;

    KinovaWaitrPybindWrapper(const std::string urdf_filename,
                             const bool display_info);

    // Destructor
    ~KinovaWaitrPybindWrapper() = default;

    // Class methods
    void set_obstacles(const nb_2d_double obstacles_inp);

    void set_contact_surface_parameters(const double mu_inp,
                                        const double R_inp,
                                        const double maxSuctionForce_inp = 0.0,
                                        const double contactForceBuffer_inp = 0.0,
                                        const double frictionForceBuffer_inp = 0.0,
                                        const double ZMPBuffer_inp = 0.0);

    void set_end_effector(const nb_1d_double contact_position,
                          const double object_mass,
                          const nb_1d_double object_com,
                          const nb_2d_double object_inertia);

    void set_ipopt_parameters(const double tol,
                              const double constr_viol_tol,
                              const double obj_scaling_factor,
                              const double max_wall_time, 
                              const int print_level,
                              const std::string mu_strategy,
                              const std::string linear_solver,
                              const bool gradient_check);

    void set_trajectory_parameters(const nb_1d_double q0_inp,
                                   const nb_1d_double q_d0_inp,
                                   const nb_1d_double q_dd0_inp,
                                   const double duration_inp);

    void set_buffer(const nb_1d_double joint_limits_buffer_inp,
                    const nb_1d_double velocity_limits_buffer_inp,
                    const nb_1d_double torque_limits_buffer_inp);

    void set_target(const nb_1d_double q_des_inp,
                    const double tplan_inp);

    nb::tuple optimize();

    nb::tuple analyze_solution();

    // Class members
    // robot model
    Model model;

    // obstacle information
    int num_obstacles = 0;
    std::vector<Vec3> boxCenters;
    std::vector<Vec3> boxOrientation;
    std::vector<Vec3> boxSize;

    // trajectory information
    ArmourTrajectoryParameters atp;
    double T = 1;
    int N = 32;
    int degree = ARMOUR_BEZIER_CURVE_DEGREE;
    VecX q_des;
    double tplan = 0;
    int tplan_n = 0;

    // contact surface parameters
    circleContactSurfaceParams csp;

    // buffer information
    VecX joint_limits_buffer;
    VecX velocity_limits_buffer;
    VecX torque_limits_buffer;
    
    SmartPtr<KinovaWaitrOptimizer> mynlp;
    SmartPtr<IpoptApplication> app;

    // results
    VecX solution;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> trajInfo;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> contactInfo;

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
