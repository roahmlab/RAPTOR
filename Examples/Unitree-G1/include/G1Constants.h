#pragma once

namespace RAPTOR {
namespace G1 {

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
    -2.5307,         // left_hip_pitch_joint
    -0.5236,         // left_hip_roll_joint
    -2.7576,         // left_hip_yaw_joint
    -0.087267,       // left_knee_joint
    -0.87267         // left_ankle_pitch_joint
    -0.2618,         // left_ankle_roll_joint
    -2.5307,         // right_hip_pitch_joint
    -2.9671,         // right_hip_roll_joint
    -2.7576,         // right_hip_yaw_joint 
    -0.087267,       // right_knee_joint   
    -0.87267,        // right_ankle_pitch_joint
    -0.2618          // right_ankle_roll_joint
};

// pulled from talos_reduced_armfixed.urdf
constexpr double JOINT_LIMITS_UPPER[NUM_JOINTS] = {
    1000,           // Px
    1000,           // Py
    1000,           // Pz
    1000,           // Rx
    1000,           // Ry
    1000,           // Rz
    2.8798,         // left_hip_pitch_joint
    2.9671,         // left_hip_roll_joint
    2.7576,         // left_hip_yaw_joint
    2.8798,         // left_knee_joint
    0.5236,         // left_ankle_pitch_joint
    0.2618,         // left_ankle_roll_joint   
    2.8798,         // right_hip_pitch_joint
    0.5236,         // right_hip_roll_joint
    2.7576,         // right_hip_yaw_joint
    2.8798,         // right_knee_joint
    0.5236,          // right_ankle_pitch_joint
    0.2618          // right_ankle_roll_joint
};

constexpr double TORQUE_LIMITS_LOWER[NUM_INDEPENDENT_JOINTS] = {
    -88,  // left_hip_pitch_joint
    -88,  // left_hip_roll_joint
    -88,  // left_hip_yaw_joint
    -139, // left_knee_joint
    -50,  // left_ankle_pitch_joint
    -50,  // left_ankle_roll_joint
    -88, // right_hip_pitch_joint
    -88, // right_hip_roll_joint
    -88, // right_hip_yaw_joint
    -139, // right_knee_joint
    -50, // right_ankle_pitch_joint
    -50  // right_ankle_roll_joint
};

constexpr double TORQUE_LIMITS_UPPER[NUM_INDEPENDENT_JOINTS] = {
    88,  // left_hip_pitch_joint
    88,  // left_hip_roll_joint
    88,  // left_hip_yaw_joint
    139, // left_knee_joint
    50,  // left_ankle_pitch_joint
    50,  // left_ankle_roll_joint
    88,  // right_hip_pitch_joint
    88,  // right_hip_roll_joint
    88,  // right_hip_yaw_joint
    139, // right_knee_joint
    50,  // right_ankle_pitch_joint
    50   // right_ankle_roll_joint
};

constexpr double MU = 0.7;
constexpr double GAMMA = 0.7;
constexpr double FOOT_WIDTH = 0.10; // (m)
constexpr double FOOT_LENGTH = 0.03; // (m)

constexpr char LEFT_FOOT_NAME[] = "left_ankle_roll_joint";
constexpr char RIGHT_FOOT_NAME[] = "right_ankle_roll_joint";

}; // namespace G1
}; // namespace RAPTOR