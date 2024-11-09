#include "RobotInfo.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

Eigen::VectorXd map_to_vector(const YAML::Node& node, int size) {
    std::vector<double> vec;
    
    try {
        vec = node.as<std::vector<double>>();
    }
    catch (const std::exception& e) {
        throw std::runtime_error("Failed to load the vector of " + node.begin()->first.as<std::string>() + ".");
    }

    if (vec.size() != size) {
        throw std::runtime_error("Size of the vector of " + node.begin()->first.as<std::string>() + " does not match the expected size.");
    }

    Eigen::VectorXd result(size);

    for (int i = 0; i < size; ++i) {
        result(i) = vec[i];
    }

    return result;
}

RobotInfo::RobotInfo(const std::string& urdf_filename,
                     const std::string& yaml_filename) {
    YAML::Node config;
    
    try {
        config = YAML::LoadFile(yaml_filename);
    }
    catch (const std::exception& e) {
        throw std::runtime_error("Failed to load the YAML file." + std::string(e.what()));
    }

    if (config.size() < 1 || !config.begin()->second.IsMap()) {
        throw std::runtime_error("Invalid YAML format. Expected a map as the first node.");
    }

    const std::string& RobotName = config.begin()->first.as<std::string>();
    const YAML::Node& RobotConfig = config[RobotName];

    try {
        pinocchio::Model model_double;
        pinocchio::urdf::buildModel(urdf_filename, model_double);
        model = model_double.cast<double>();
    }
    catch (const std::exception& e) {
        throw std::runtime_error("Failed to load the URDF file: " + std::string(e.what()));
    }

    model.name = RobotName;

    num_motors = RobotConfig["num_motors"].as<int>();
    num_joints = RobotConfig["num_joints"].as<int>();

    if (num_motors != NUM_FACTORS) {
        std::cerr << "For now we only support 7 variables and the number is hardcoded unfortunately.\n";
        std::cerr << "The degree of the robot has to be 7.\n";
        std::cerr << "If not, we need to change the following macros to support a different number of variables.\n";
        throw std::runtime_error("Number of motors in the YAML file does not match the number of motors in the code.");
    }

    if (model.nv != num_joints) {
        throw std::runtime_error("Number of joints in the URDF file does not match the number of joints in the YAML file.");
    }

    if (num_joints < num_motors) {
        throw std::runtime_error("Number of joints is less than the number of motors.");
    }

    Vec3 gravity = map_to_vector(RobotConfig["gravity"], 3);
    model.gravity.linear() = gravity;

    model.friction.setZero();
    model.friction.head(num_motors) = map_to_vector(RobotConfig["friction"], num_motors);
    model.damping.setZero();
    model.damping.head(num_motors) = map_to_vector(RobotConfig["damping"], num_motors);
    model.armature.setZero();
    model.armature.head(num_motors) = map_to_vector(RobotConfig["transmissionInertia"], num_motors);
    position_limits_lb = map_to_vector(RobotConfig["position_limits_lb"], num_motors);
    position_limits_ub = map_to_vector(RobotConfig["position_limits_ub"], num_motors);
    velocity_limits = map_to_vector(RobotConfig["velocity_limits"], num_motors);
    torque_limits = map_to_vector(RobotConfig["torque_limits"], num_motors);
    mass_uncertainty = map_to_vector(RobotConfig["mass_uncertainty"], num_joints);
    inertia_uncertainty = map_to_vector(RobotConfig["inertia_uncertainty"], num_joints);

    try {
        const YAML::Node& ultimate_bound_node = RobotConfig["ultimate_bound"];
        ultimate_bound_info.alpha = ultimate_bound_node["alpha"].as<double>();
        ultimate_bound_info.V_m = ultimate_bound_node["V_m"].as<double>();
        ultimate_bound_info.M_max = ultimate_bound_node["M_max"].as<double>();
        ultimate_bound_info.M_min = ultimate_bound_node["M_min"].as<double>();
        ultimate_bound_info.Kr = ultimate_bound_node["Kr"].as<double>();

        ultimate_bound_info.eps = std::sqrt(2 * ultimate_bound_info.V_m / ultimate_bound_info.M_min);
        ultimate_bound_info.qe = ultimate_bound_info.eps / ultimate_bound_info.Kr;
        ultimate_bound_info.qde = 2 * ultimate_bound_info.eps;
        ultimate_bound_info.qdae = ultimate_bound_info.eps;
        ultimate_bound_info.qddae = 2 * ultimate_bound_info.Kr * ultimate_bound_info.eps;
    }
    catch (const std::exception& e) {
        throw std::runtime_error("Failed to load the ultimate bound information.");
    }

    num_spheres = 0;
    for (const auto& entry : RobotConfig["collision_spheres"]) {
        std::string link_name = entry.first.as<std::string>();
        const YAML::Node& spheres = entry.second;

        std::vector<std::pair<Vec3, double>> spheres_info;

        for (const auto& sphere : spheres) {
            const YAML::Node& offset_node = sphere["offset"];

            Vec3 offset;
            offset << offset_node[0].as<double>(), offset_node[1].as<double>(), offset_node[2].as<double>();
            double radius = sphere["radius"].as<double>();

            spheres_info.emplace_back(offset, radius);
        }

        collision_spheres[link_name] = spheres_info;
        num_spheres += spheres_info.size();
    }
}

void RobotInfo::print() const {
    std::cout << "Model: " << model.name << std::endl;
    std::cout << "Number of joints: " << num_joints << std::endl;
    std::cout << "Number of motors: " << num_motors << std::endl;
    std::cout << "Gravity:\n    " << model.gravity.linear().transpose() << std::endl;
    std::cout << "Friction:\n    " << model.friction.transpose() << std::endl;
    std::cout << "Damping:\n    " << model.damping.transpose() << std::endl;
    std::cout << "Rotor Inertia:\n    " << model.armature.transpose() << std::endl;
    std::cout << "Position Limits Lower Bound:\n    " << position_limits_lb.transpose() << std::endl;
    std::cout << "Position Limits Upper Bound:\n    " << position_limits_ub.transpose() << std::endl;
    std::cout << "Velocity Limits:\n    " << velocity_limits.transpose() << std::endl;
    std::cout << "Torque Limits:\n    " << torque_limits.transpose() << std::endl;
    std::cout << "Ultimate Bound Information:" << std::endl;
    std::cout << "    Alpha: " << ultimate_bound_info.alpha << std::endl;
    std::cout << "    V_m: " << ultimate_bound_info.V_m << std::endl;
    std::cout << "    M_max: " << ultimate_bound_info.M_max << std::endl;
    std::cout << "    M_min: " << ultimate_bound_info.M_min << std::endl;
    std::cout << "    Kr: " << ultimate_bound_info.Kr << std::endl;
    std::cout << "Collision Spheres:" << std::endl;
    for (const auto& entry : collision_spheres) {
        std::cout << "    " << entry.first << std::endl;
        for (const auto& sphere : entry.second) {
            std::cout << "        Offset: " << sphere.first.transpose() << std::endl;
            std::cout << "        Radius: " << sphere.second << std::endl;
        }
    }
}

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR