#ifndef ROBOTINFO_H
#define ROBOTINFO_H

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include "PZsparse.h"

#include <unordered_map>
#include <yaml-cpp/yaml.h>

namespace RAPTOR {
namespace Armour {
    
struct ultimate_bound {
    // controller information
    float alpha = 0;
    float V_m = 0;
    float M_max = 0;
    float M_min = 0;
    float Kr = 0;

    // ultimate bound information
    float eps = 0;
    float qe = 0;
    float qde = 0;
    float qdae = 0;
    float qddae = 0;
};

Eigen::VectorXf map_to_vector(const YAML::Node& node, int size);

class RobotInfo {
public:
    using Vec3 = Eigen::Vector3f;
    using VecX = Eigen::VectorXf;

    int num_motors = 0;
    int num_joints = 0;

    pinocchio::ModelTpl<float> model;

    VecX mass_uncertainty;
    VecX inertia_uncertainty;

    VecX position_limits_lb;
    VecX position_limits_ub;
    VecX velocity_limits;
    VecX torque_limits;

    ultimate_bound ultimate_bound_info;

    int num_spheres = 0;
    std::unordered_map<std::string, std::vector<std::pair<Vec3, float>>> collision_spheres;

    RobotInfo() = default;

    RobotInfo(const std::string& urdf_filename,
              const std::string& yaml_filename);

    ~RobotInfo() = default;

    void print() const;
};

}; // namespace Armour
}; // namespace RAPTOR

#endif // ROBOTINFO_H