#ifndef DUAL_ARMOUR_PYBIND_WRAPPER_H
#define DUAL_ARMOUR_PYBIND_WRAPPER_H

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include "DualArmourOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

namespace nb = nanobind;

class DualArmourPybindWrapper {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    using nb_1d_double = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
    using nb_2d_double = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

    // Constructor
    DualArmourPybindWrapper() = default;

    DualArmourPybindWrapper(const std::string urdf_filename1,
                            const std::string config_filename1,
                            const std::string urdf_filename2,
                            const std::string config_filename2,
                            const bool display_info = false);

    // Destructor
    ~DualArmourPybindWrapper() = default;

    // Class methods
    void set_obstacles(const nb_2d_double obstacles_inp);

    void set_endeffector_inertial_parameters_robot1(const double object_mass,
                                                    const nb_1d_double object_com,
                                                    const nb_1d_double object_inertia);

    void set_endeffector_inertial_parameters_robot2(const double object_mass,
                                                    const nb_1d_double object_com,
                                                    const nb_1d_double object_inertia);

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
                                   const nb_1d_double k_center_inp,
                                   const nb_1d_double k_range_inp,
                                   const double duration_inp,
                                   const nb_1d_double q_des_inp,
                                   const double t_plan_inp);

    nb::tuple optimize();

    nb::tuple analyze_solution();

    // Class members
    std::shared_ptr<RobotInfo> robotInfoPtr1_ = nullptr;
    std::shared_ptr<RobotInfo> robotInfoPtr2_ = nullptr;

    std::shared_ptr<BezierCurveInterval> trajPtr1_ = nullptr;
    std::shared_ptr<BezierCurveInterval> trajPtr2_ = nullptr;

    std::shared_ptr<PZDynamics> dynPtr1_ = nullptr;
    std::shared_ptr<PZDynamics> dynPtr2_ = nullptr;

    // obstacle information
    int num_obstacles = 0;
    std::vector<Vec3> boxCenters;
    std::vector<Vec3> boxOrientation;
    std::vector<Vec3> boxSize;

    // trajectory parameters
    VecX q0_robot1;
    VecX q_d0_robot1;
    VecX q_dd0_robot1;
    VecX k_center_robot1;
    VecX k_range_robot1;
    VecX q0_robot2;
    VecX q_d0_robot2;
    VecX q_dd0_robot2;
    VecX k_center_robot2;
    VecX k_range_robot2;
    double duration = 3.0;
    VecX q_des;
    double t_plan = 0.0;
    
    SmartPtr<DualArmourOptimizer> mynlp;
    SmartPtr<IpoptApplication> app;

    // detailed information of the solution trajectory and the corresponding reachable sets
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> trajInfo_robot1;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> trajInfo_robot2;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> spheres_x;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> spheres_y;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> spheres_z;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> spheres_radius;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> torque_center;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> torque_radius;
    // Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> separation_force_center;
    // Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> separation_force_radius;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> friction_cone_center;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> friction_cone_radius;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> zmp_center;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> zmp_radius;

    // Flags to check if the parameters are set
    bool set_obstacles_check = false;
    bool set_ipopt_parameters_check = false;
    bool set_trajectory_parameters_check = false;
    bool has_optimized = false;
};

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR

#endif // DUAL_ARMOUR_PYBIND_WRAPPER_H
