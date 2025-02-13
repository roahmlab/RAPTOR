#ifndef KINOVA_IK_MOTION_PYBIND_WRAPPER_H

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include "KinovaIKSolver.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

namespace RAPTOR {
namespace Kinova {

namespace nb = nanobind;

class KinovaIKMotionPybindWrapper {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using Mat3 = Eigen::Matrix3d;
    using MatX = Eigen::MatrixXd;

    using nb_1d_double = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
    using nb_2d_double = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

    // Constructor
    KinovaIKMotionPybindWrapper() = default;

    KinovaIKMotionPybindWrapper(const std::string urdf_filename,
                                const bool display_info);

    // Destructor
    ~KinovaIKMotionPybindWrapper() = default;

    // Class methods
    void set_desired_endeffector_transforms(const nb_2d_double& desired_endeffector_transforms_inp);

    void set_obstacles(const nb_2d_double& obstacles_inp,
                       const double collision_buffer_inp);

    void set_ipopt_parameters(const double tol,
                              const double constr_viol_tol,
                              const double obj_scaling_factor,
                              const double max_wall_time, 
                              const int print_level,
                              const std::string mu_strategy,
                              const std::string linear_solver,
                              const bool gradient_check);

    nb::tuple solve(const nb_1d_double& initial_guess);

    // Class members
    // robot model
    Model model;

    // end effector transform
    Transform endT;

    // obstacle information
    int num_obstacles = 0;
    std::vector<Vec3> boxCenters;
    std::vector<Vec3> boxOrientation;
    std::vector<Vec3> boxSize;
    double collision_buffer = 0;
    
    SmartPtr<KinovaIKSolver> mynlp;
    SmartPtr<IpoptApplication> app;

    std::vector<Transform> desiredTransforms;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> solutions;

    // Flags to check if the parameters are set
    bool set_transform_check = false;
    bool set_ipopt_parameters_check = false;
    bool has_optimized = false;
};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // KINOVA_IK_MOTION_PYBIND_WRAPPER_H