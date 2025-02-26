#pragma once

namespace RAPTOR {
namespace Digit {

// constants related to digit-v3
constexpr int NUM_JOINTS = 36;
constexpr int NUM_DEPENDENT_JOINTS = 24;
constexpr int NUM_INDEPENDENT_JOINTS = 12;

// pulled out from digit-v3.xml
// This is in degree!!!
constexpr double JOINT_LIMITS_LOWER[NUM_JOINTS] = {
    -1000,    // Px
    -1000,    // Py
    -1000,    // Pz
    -1000,    // Rx
    -1000,    // Ry
    -1000,    // Rz
    -60,      // left-hip-roll
    -40,      // left-hip-yaw
    -60,      // left-hip-pitch
    -180,     // left-achilles-rod
    -10,      // left-ach2
    -16.5480 - 15, // -80,      // left-knee
    -50.3,    // left-tarsus
    -26.8765 - 15, // -46.2755, // left-toe-A
    -180,     // left-toe-A-rod
    -10,      // left-A2
    -1.1295 - 15, // -45.8918, // left-toe-B
    -180,     // left-toe-B-rod
    -10,      // left-B2
    -44,      // left-toe-pitch
    -37,      // left-toe-roll
    -60,      // right-hip-roll
    -40,      // right-hip-yaw
    -90,      // right-hip-pitch
    -180,     // right-achilles-rod
    -10,      // right-ach2
    -58.4,    // right-knee
    -71.6,    // right-tarsus
    -18.2531 - 15, // -44.9815, // right-toe-A
    -180,     // right-toe-A-rod
    -10,      // right-A2
    -26.2724 - 15, // -45.5476, // right-toe-B
    -180,     // right-toe-B-rod
    -10,      // right-B2
    -34,      // right-toe-pitch
    -33       // right-toe-roll
};

// pulled out from digit-v3.xml
// This is in degree!!!
constexpr double JOINT_LIMITS_UPPER[NUM_JOINTS] = {
    1000,    // Px
    1000,    // Py
    1000,    // Pz
    1000,    // Rx
    1000,    // Ry
    1000,    // Rz
    60,      // left-hip-roll
    40,      // left-hip-yaw
    90,      // left-hip-pitch
    180,     // left-achilles-rod
    10,      // left-ach2
    58.4,    // left-knee
    71.6,    // left-tarsus
    0.4032 + 15, // 44.9815, // left-toe-A
    180,     // left-toe-A-rod
    10,      // left-A2
    26.2183 + 15, // 45.5476, // left-toe-B
    180,     // left-toe-B-rod
    10,      // left-B2
    34,      // left-toe-pitch
    33,      // left-toe-roll
    60,      // right-hip-roll
    40,      // right-hip-yaw
    60,      // right-hip-pitch
    180,     // right-achilles-rod
    10,      // right-ach2
    30.7455 + 15, // 80,      // right-knee
    50.3,    // right-tarsus
    27.0323 + 15, // 46.2755, // right-toe-A
    180,     // right-toe-A-rod
    10,      // right-A2
    14.6109 + 15, // 45.8918, // right-toe-B
    180,     // right-toe-B-rod
    10,      // right-B2
    44,      // right-toe-pitch
    37       // right-toe-roll
};

// pulled out from digit-v3.xml
constexpr double TORQUE_LIMITS_LOWER[NUM_INDEPENDENT_JOINTS] = {
    -113.000,   // left-hip-roll
    -79.1765,   // left-hip-yaw
    -216.928,   // left-hip-pitch
    -231.317,   // left-knee
    -41.9759,   // left-toe-A
    -41.9759,   // left-toe-B
    -126.682,   // right-hip-roll
    -79.1765,   // right-hip-yaw
    -216.928,   // right-hip-pitch
    -231.317,   // right-knee
    -41.9759,   // right-toe-A
    -41.9759,   // right-toe-B
};

// pulled out from digit-v3.xml
constexpr double TORQUE_LIMITS_UPPER[NUM_INDEPENDENT_JOINTS] = {
    113.000,   // left-hip-roll
    79.1765,   // left-hip-yaw
    216.928,   // left-hip-pitch
    231.317,   // left-knee
    41.9759,   // left-toe-A
    41.9759,   // left-toe-B
    126.682,   // right-hip-roll
    79.1765,   // right-hip-yaw
    216.928,   // right-hip-pitch
    231.317,   // right-knee
    41.9759,   // right-toe-A
    41.9759,   // right-toe-B
};

// pulled out from digit-v3.xml
constexpr double MU = 0.7;
constexpr double GAMMA = 0.7;
constexpr double FOOT_WIDTH = 0.04; // (m)
constexpr double FOOT_LENGTH = 0.1175; // (m)

// pulled out from digit-v3.xml
constexpr double GRAVITY = -9.806; // m/s^2

constexpr char LEFT_FOOT_NAME[] = "left_toe_roll";
constexpr char RIGHT_FOOT_NAME[] = "right_toe_roll";

}; // namespace Digit
}; // namespace RAPTOR