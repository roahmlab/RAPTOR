#pragma once

namespace RAPTOR {
namespace Talos {

// constants related to talos
constexpr int NUM_JOINTS = 18;
constexpr int NUM_DEPENDENT_JOINTS = 6;
constexpr int NUM_INDEPENDENT_JOINTS = 12;

// pulled from talos_reduced_armfixed.urdf
constexpr double JOINT_LIMITS_LOWER[NUM_JOINTS] = {
    -1000,           // Px
    -1000,           // Py
    -1000,           // Pz
    -1000,           // Rx
    -1000,           // Ry
    -1000,           // Rz
    -0.349065850399, // leg_left_1_joint
    -0.5236,         // leg_left_2_joint
    -2.095,          // leg_left_3_joint
    0,               // leg_left_4_joint
    -1.309,          // leg_left_5_joint
    -0.5236,         // leg_left_6_joint
    -1.57079632679,  // leg_right_1_joint
    -0.5236,         // leg_right_2_joint
    -2.095,          // leg_right_3_joint 
    0,               // leg_right_4_joint   
    -1.309,          // leg_right_5_joint
    -0.5236          // leg_right_6_joint
};

// pulled from talos_reduced_armfixed.urdf
constexpr double JOINT_LIMITS_UPPER[NUM_JOINTS] = {
    1000,           // Px
    1000,           // Py
    1000,           // Pz
    1000,           // Rx
    1000,           // Ry
    1000,           // Rz
    1.57079632679,  // leg_left_1_joint
    0.5236,         // leg_left_2_joint
    0.7,            // leg_left_3_joint
    2.618,          // leg_left_4_joint
    0.768,          // leg_left_5_joint
    0.5236,         // leg_left_6_joint   
    0.349065850399, // leg_right_1_joint
    0.5236,         // leg_right_2_joint
    0.7,            // leg_right_3_joint
    2.618,          // leg_right_4_joint
    0.768,          // leg_right_5_joint
    0.5236          // leg_right_6_joint
};

constexpr double TORQUE_LIMITS_LOWER[NUM_INDEPENDENT_JOINTS] = {
    -100, // leg_left_1_joint
    -160, // leg_left_2_joint
    -160, // leg_left_3_joint
    -300, // leg_left_4_joint
    -160, // leg_left_5_joint
    -100, // leg_left_6_joint
    -100, // leg_right_1_joint
    -160, // leg_right_2_joint
    -160, // leg_right_3_joint
    -300, // leg_right_4_joint
    -160, // leg_right_5_joint
    -100  // leg_right_6_joint
};

constexpr double TORQUE_LIMITS_UPPER[NUM_INDEPENDENT_JOINTS] = {
    100, // leg_left_1_joint
    160, // leg_left_2_joint
    160, // leg_left_3_joint
    300, // leg_left_4_joint
    160, // leg_left_5_joint
    100, // leg_left_6_joint
    100, // leg_right_1_joint
    160, // leg_right_2_joint
    160, // leg_right_3_joint
    300, // leg_right_4_joint
    160, // leg_right_5_joint
    100  // leg_right_6_joint
};

constexpr double MU = 0.7;
constexpr double GAMMA = 0.7;
constexpr double FOOT_WIDTH = 0.1344; // (m)
constexpr double FOOT_LENGTH = 0.2208; // (m)

constexpr char LEFT_FOOT_NAME[] = "leg_left_6_joint";
constexpr char RIGHT_FOOT_NAME[] = "leg_right_6_joint";

}; // namespace Talos
}; // namespace RAPTOR