#ifndef KINOVA_PYBIND_WRAPPER_H
#define KINOVA_PYBIND_WRAPPER_H

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include "SafePayloadExcitingTrajectoryGenerator.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

namespace nb = nanobind;

class SafePayloadExcitingPybindWrapper {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using Vec10 = Eigen::Vector<double, 10>;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    using nb_1d_double = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
    using nb_2d_double = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

    // Constructor
    SafePayloadExcitingPybindWrapper() = default;

    SafePayloadExcitingPybindWrapper(const std::string urdf_filename,
                                     const std::string config_filename,
                                     const bool display_info = false);

    // Destructor
    ~SafePayloadExcitingPybindWrapper() = default;

    // Class methods
    void set_obstacles(const nb_2d_double obstacles_inp);

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
                                   const double duration_inp);

    void set_endeffector_inertial_parameters(const nb_1d_double inertial_parameters,
                                             const nb_1d_double inertial_parameters_lb,
                                             const nb_1d_double inertial_parameters_ub);

    nb::tuple optimize();

    nb::tuple analyze_solution();

    // Class members
    std::shared_ptr<RobotInfo> robotInfoPtr_ = nullptr;

    std::shared_ptr<BezierCurveInterval> trajPtr_ = nullptr;
    std::shared_ptr<PZDynamics> dynPtr_ = nullptr;

    // obstacle information
    int num_obstacles = 0;
    std::vector<Vec3> boxCenters;
    std::vector<Vec3> boxOrientation;
    std::vector<Vec3> boxSize;

    // trajectory parameters
    VecX q0;
    VecX q_d0;
    VecX q_dd0;
    VecX k_center;
    VecX k_range;
    double duration;
    
    SmartPtr<SafePayloadExcitingTrajectoryGenerator> mynlp;
    SmartPtr<IpoptApplication> app;

    // detailed information of the solution trajectory and the corresponding reachable sets
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> trajInfo;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> spheres_x;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> spheres_y;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> spheres_z;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> spheres_radius;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> torque_center;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> torque_radius;

    // Flags to check if the parameters are set
    bool set_obstacles_check = false;
    bool set_ipopt_parameters_check = false;
    bool set_trajectory_parameters_check = false;
    bool has_optimized = false;
};

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR

#endif // KINOVA_PYBIND_WRAPPER_H
