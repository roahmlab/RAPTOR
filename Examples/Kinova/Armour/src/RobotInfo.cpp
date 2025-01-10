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
    mass_uncertainty = map_to_vector(RobotConfig["mass_uncertainty"], num_joints);
    com_uncertainty = map_to_vector(RobotConfig["com_uncertainty"], num_joints);
    inertia_uncertainty = map_to_vector(RobotConfig["inertia_uncertainty"], num_joints);

    if (RobotConfig["end_effector_mass_lb"].IsDefined() &&
        RobotConfig["end_effector_mass_ub"].IsDefined() &&
        RobotConfig["end_effector_com_lb"].IsDefined() &&
        RobotConfig["end_effector_com_ub"].IsDefined() &&
        RobotConfig["end_effector_inertia_lb"].IsDefined() &&
        RobotConfig["end_effector_inertia_ub"].IsDefined()) {
        if_end_effector_info_exist = true;
        end_effector_mass_lb = RobotConfig["end_effector_mass_lb"].as<double>();
        end_effector_mass_ub = RobotConfig["end_effector_mass_ub"].as<double>();
        end_effector_com_lb = map_to_vector(RobotConfig["end_effector_com_lb"], 3);
        end_effector_com_ub = map_to_vector(RobotConfig["end_effector_com_ub"], 3);
        end_effector_inertia_lb = map_to_vector(RobotConfig["end_effector_inertia_lb"], 6);
        end_effector_inertia_ub = map_to_vector(RobotConfig["end_effector_inertia_ub"], 6);

        // consistency check
        const double mass = model.inertias[model.nv].mass();
        if (end_effector_mass_lb > mass ||
            end_effector_mass_ub < mass) {
            std::cerr << mass << " [" << end_effector_mass_lb << ", " << end_effector_mass_ub << "]\n";
            throw std::runtime_error("End effector mass bounds are not within the range of the original mass.");
        }
        const Vec3& com = model.inertias[model.nv].lever();
        if ((end_effector_com_lb - com).minCoeff() > 0 ||
            (end_effector_com_ub - com).maxCoeff() < 0) {
            std::cerr << com.transpose() << "\n[" << end_effector_com_lb.transpose() << ",\n " << end_effector_com_ub.transpose() << "]\n";
            throw std::runtime_error("End effector com bounds are not within the range of the original com.");
        }
        const Vec6& inertia = model.inertias[model.nv].inertia().data();
        if ((end_effector_inertia_lb - inertia).minCoeff() > 0 ||
            (end_effector_inertia_ub - inertia).maxCoeff() < 0) {
            std::cerr << inertia.transpose() << "\n[" << end_effector_inertia_lb.transpose() << ",\n " << end_effector_inertia_ub.transpose() << "]\n";
            throw std::runtime_error("End effector inertia bounds are not within the range of the original inertia.");
        }
    }

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

    if (num_joints > num_motors) { // solving waitr problem, need information of the contact surface
        try {
            suction_force = RobotConfig["suction_force"].as<double>();
            mu = RobotConfig["mu"].as<double>();
            contact_surface_radius = RobotConfig["contact_surface_radius"].as<double>();
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Failed to load the contact surface information.");
        }
    }

    num_spheres = 0;
    num_capsules = 0;
    sphere_radii.clear();
    for (const auto& entry : RobotConfig["collision_spheres"]) {
        std::string link_name = entry.first.as<std::string>();
        const YAML::Node& spheres = entry.second;
        
        if (model.existJointName(link_name)) {
            for (const auto& sphere : spheres) {
                const YAML::Node& offset_node = sphere["offset"];
                const YAML::Node& radius_node = sphere["radius"];
                Vec3 offset;
                offset << offset_node[0].as<double>(), 
                          offset_node[1].as<double>(), 
                          offset_node[2].as<double>();

                sphere_radii.push_back(radius_node.as<double>());

                pinocchio::SE3 sphere_placement;
                sphere_placement.setIdentity();
                sphere_placement.translation() = offset;

                model.addFrame(
                    pinocchio::Frame(
                        "collision-" + std::to_string(num_spheres),
                        model.getJointId(link_name),
                        0,
                        sphere_placement,
                        pinocchio::OP_FRAME)
                    );

                num_spheres++;
            }
        }
        else {
            throw std::runtime_error("Link " + link_name + " does not exist in the URDF file.");
        }
    }

    // Import tapered capsules
    for (const auto& entry : RobotConfig["tapered_capsules"]){
        std::string link_name = entry.first.as<std::string>();
        const YAML::Node& spheres = entry.second;
        for (const auto& sphere : spheres) {
            const size_t sphere_1 = sphere["sphere_1"].as<size_t>();
            const size_t sphere_2 = sphere["sphere_2"].as<size_t>();
            tc_begin_and_end.push_back(std::make_pair(sphere_1, sphere_2));
            num_capsules++;
        }
    }

    // Generate list of self-collision checks
    num_self_collisions = 0;
    for (int arm_1_index = 0; arm_1_index < num_capsules - 2; arm_1_index++){
        for (int arm_2_index = arm_1_index + 2; arm_2_index < num_capsules; arm_2_index++){
            self_collision_checks.push_back(std::make_pair(arm_1_index, arm_2_index));
            num_self_collisions++;
        }
    }
}

void RobotInfo::change_endeffector_inertial_parameters(const double mass,
                                                       const Vec3& com,
                                                       const Vec6& inertia,
                                                       const double mass_eps,
                                                       const double com_eps,
                                                       const double inertia_eps) {
    Mat3 inertia_matrix;
    inertia_matrix << inertia(0), inertia(1), inertia(3),
                      inertia(1), inertia(2), inertia(4),
                      inertia(3), inertia(4), inertia(5);
    model.inertias[model.nv - 1] = pinocchio::Inertia(mass, com, inertia_matrix);

    if (mass_eps < 0 ||
        com_eps < 0 ||
        inertia_eps < 0) {
        throw std::invalid_argument("Uncertainty values cannot be negative.");
    }

    mass_uncertainty(model.nv - 1) = mass_eps;
    com_uncertainty(model.nv - 1)= com_eps;
    inertia_uncertainty(model.nv - 1) = inertia_eps;
}

void RobotInfo::change_endeffector_inertial_parameters(const Vec10& inertial_parameters,
                                                       const Vec10& inertial_parameters_lb,
                                                       const Vec10& inertial_parameters_ub) {
    model.inertias[model.nv] = pinocchio::Inertia::FromDynamicParameters(inertial_parameters);

    if_end_effector_info_exist = true;

    end_effector_mass_lb = inertial_parameters_lb(0);
    end_effector_mass_ub = inertial_parameters_ub(0);

    end_effector_com_lb = inertial_parameters_lb.segment<3>(1);
    end_effector_com_ub = inertial_parameters_ub.segment<3>(1);

    end_effector_inertia_lb = inertial_parameters_lb.segment<6>(4);
    end_effector_inertia_ub = inertial_parameters_ub.segment<6>(4);

    // consistency check
    const double mass = model.inertias[model.nv].mass();
    if (end_effector_mass_lb > mass ||
        end_effector_mass_ub < mass) {
        std::cerr << mass << " [" << end_effector_mass_lb << ", " << end_effector_mass_ub << "]\n";
        throw std::runtime_error("End effector mass bounds are not within the range of the original mass.");
    }
    const Vec3& com = model.inertias[model.nv].lever();
    if ((end_effector_com_lb - com).minCoeff() > 0 ||
        (end_effector_com_ub - com).maxCoeff() < 0) {
        std::cerr << com.transpose() << "\n[" << end_effector_com_lb.transpose() << ",\n " << end_effector_com_ub.transpose() << "]\n";
        throw std::runtime_error("End effector com bounds are not within the range of the original com.");
    }
    const Vec6& inertia = model.inertias[model.nv].inertia().data();
    if ((end_effector_inertia_lb - inertia).minCoeff() > 0 ||
        (end_effector_inertia_ub - inertia).maxCoeff() < 0) {
        std::cerr << inertia.transpose() << "\n[" << end_effector_inertia_lb.transpose() << ",\n " << end_effector_inertia_ub.transpose() << "]\n";
        throw std::runtime_error("End effector inertia bounds are not within the range of the original inertia.");
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
    std::cout << "Ultimate Bound Information:" << std::endl;
    std::cout << "    Alpha: " << ultimate_bound_info.alpha << std::endl;
    std::cout << "    V_m: " << ultimate_bound_info.V_m << std::endl;
    std::cout << "    M_max: " << ultimate_bound_info.M_max << std::endl;
    std::cout << "    M_min: " << ultimate_bound_info.M_min << std::endl;
    std::cout << "    Kr: " << ultimate_bound_info.Kr << std::endl;
}

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR