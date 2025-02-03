#ifndef KINOVA_LONGER_HORIZON_PYBIND_WRAPPER_H
#define KINOVA_LONGER_HORIZON_PYBIND_WRAPPER_H

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include "KinovaLongerHorizonOptimizer.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

namespace RAPTOR {
namespace Kinova {

namespace nb = nanobind;

class KinovaLongerHorizonPybindWrapper {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    using nb_1d_double = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
    using nb_2d_double = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

    // Constructor
    KinovaLongerHorizonPybindWrapper() = default;

    KinovaLongerHorizonPybindWrapper(const std::string urdf_filename,
                                     const bool display_info);

    // Destructor
    ~KinovaLongerHorizonPybindWrapper() = default;

    // Class methods
    void set_obstacles(const nb_2d_double obstacles_inp,
                       const double collision_buffer_inp);

    void set_ipopt_parameters(const double tol,
                              const double constr_viol_tol,
                              const double obj_scaling_factor,
                              const double max_wall_time, 
                              const int print_level,
                              const std::string mu_strategy,
                              const std::string linear_solver,
                              const bool gradient_check);

    void set_trajectory_parameters(const nb_1d_double q0_inp,
                                   const nb_1d_double qT_inp, 
                                   const double duration_inp,
                                   const int degree_inp);

    void set_buffer(const nb_1d_double joint_limits_buffer_inp,
                    const nb_1d_double velocity_limits_buffer_inp,
                    const nb_1d_double torque_limits_buffer_inp);

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
    double collision_buffer = 0;

    // trajectory information
    double T = 5.0;
    int N = 64;
    int degree = 3;
    VecX q0;
    VecX qT;

    // buffer information
    VecX joint_limits_buffer;
    VecX velocity_limits_buffer;
    VecX torque_limits_buffer;
    
    SmartPtr<KinovaLongerHorizonOptimizer> mynlp;
    SmartPtr<IpoptApplication> app;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> trajInfo; // trajectory information
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> spheres_x;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> spheres_y;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> spheres_z;

    // Flags to check if the parameters are set
    bool set_obstacles_check = false;
    bool set_ipopt_parameters_check = false;
    bool set_trajectory_parameters_check = false;
    bool set_buffer_check = false;
    bool has_optimized = false;
};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // KINOVA_LONGER_HORIZON_PYBIND_WRAPPER_H
