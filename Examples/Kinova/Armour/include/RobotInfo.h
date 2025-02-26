#ifndef ROBOTINFO_H
#define ROBOTINFO_H

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include "KinovaConstants.h"

#include "PZSparse.h"

#include <unordered_map>
#include <yaml-cpp/yaml.h>

namespace RAPTOR {
namespace Kinova {
namespace Armour {
    
struct ultimate_bound {
    // controller information
    double alpha = 0;
    double V_m = 0;
    double M_max = 0;
    double M_min = 0;
    double Kr = 0;

    // ultimate bound information
    double eps = 0;
    double qe = 0;
    double qde = 0;
    double qdae = 0;
    double qddae = 0;
};

Eigen::VectorXd map_to_vector(const YAML::Node& node, int size);

class RobotInfo {
public:
    using Vec3 = Eigen::Vector3d;
    using Vec6 = Eigen::Vector<double, 6>;
    using Vec10 = Eigen::Vector<double, 10>;
    using VecX = Eigen::VectorXd;
    using Mat3 = Eigen::Matrix3d;

    int num_motors = 0;
    int num_joints = 0;

    pinocchio::Model model;

    // model uncertainty
    VecX mass_uncertainty;
    VecX com_uncertainty;
    VecX inertia_uncertainty;

    // the following are specifically for the end effector
    // if defined, these will overwrite the previous model uncertainty 
    // (for example, the last element of the mass_uncertainty vector is the end effector mass uncertainty,
    // but will be ignored if end_effector_mass_lb and end_effector_mass_ub are defined)
    bool if_end_effector_info_exist = false;
    double end_effector_mass_lb = 0;
    double end_effector_mass_ub = 0;
    Vec3 end_effector_com_lb;
    Vec3 end_effector_com_ub;
    Vec6 end_effector_inertia_lb;
    Vec6 end_effector_inertia_ub;

    ultimate_bound ultimate_bound_info;

    // contact surface infomation
    double suction_force = 0; // additional force applied to the end effector when using suction cup
    double mu = 0; // friction coefficient of the contact surface
    double contact_surface_radius = 0; // radius of the contact surface (radius of the suction cup)

    // collision spheres info
    int num_spheres = 0;
    std::vector<double> sphere_radii;

    // bimanual constraints info
    int num_self_collisions = 0;
    int num_capsules = 0;
    std::vector<std::pair<size_t, size_t>> tc_begin_and_end; // index of the two spheres that consist a tapered capsule
    std::vector<std::pair<size_t, size_t>> self_collision_checks;

    RobotInfo() = default;

    RobotInfo(const std::string& urdf_filename,
              const std::string& yaml_filename);

    ~RobotInfo() = default;

    void change_endeffector_inertial_parameters(const double mass,
                                                const Vec3& com,
                                                const Vec6& inertia,
                                                const double mass_eps,
                                                const double com_eps,
                                                const double inertia_eps);

    void change_endeffector_inertial_parameters(const Vec10& inertial_parameters,
                                                const Vec10& inertial_parameters_lb,
                                                const Vec10& inertial_parameters_ub);

    void print() const;
};

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR

#endif // ROBOTINFO_H