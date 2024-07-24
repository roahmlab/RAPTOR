#pragma once

namespace RAPTOR {
namespace DigitModified {

// constants related to digit-v3-modified
constexpr int NUM_JOINTS = 20;
constexpr int NUM_DEPENDENT_JOINTS = 6;
constexpr int NUM_INDEPENDENT_JOINTS = 14;

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
    -80,      // left-knee
    -50.3,    // left-tarsus
    -44,      // left-toe-pitch
    -37,      // left-toe-roll
    -60,      // right-hip-roll
    -40,      // right-hip-yaw
    -90,      // right-hip-pitch
    -58.4,    // right-knee
    -71.6,    // right-tarsus
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
    58.4,    // left-knee
    71.6,    // left-tarsus
    34,      // left-toe-pitch
    33,      // left-toe-roll
    60,      // right-hip-roll
    40,      // right-hip-yaw
    60,      // right-hip-pitch
    80,      // right-knee
    50.3,    // right-tarsus
    44,      // right-toe-pitch
    37       // right-toe-roll
};

// // pulled out from digit-v3.xml
// constexpr double TORQUE_LIMITS_LOWER[NUM_INDEPENDENT_JOINTS] = {
//     -126.682,   // left-hip-roll
//     -79.1765,   // left-hip-yaw
//     -216.928,   // left-hip-pitch
//     -231.317,   // left-knee
//     -41.9759,   // left-toe-A
//     -41.9759,   // left-toe-B
//     -126.682,   // right-hip-roll
//     -79.1765,   // right-hip-yaw
//     -216.928,   // right-hip-pitch
//     -231.317,   // right-knee
//     -41.9759,   // right-toe-A
//     -41.9759,   // right-toe-B
// };

// // pulled out from digit-v3.xml
// constexpr double TORQUE_LIMITS_UPPER[NUM_INDEPENDENT_JOINTS] = {
//     126.682,   // left-hip-roll
//     79.1765,   // left-hip-yaw
//     216.928,   // left-hip-pitch
//     231.317,   // left-knee
//     41.9759,   // left-toe-A
//     41.9759,   // left-toe-B
//     126.682,   // right-hip-roll
//     79.1765,   // right-hip-yaw
//     216.928,   // right-hip-pitch
//     231.317,   // right-knee
//     41.9759,   // right-toe-A
//     41.9759,   // right-toe-B
// };

// pulled out from digit-v3.xml
constexpr double MU = 0.7;
constexpr double GAMMA = 0.7;
constexpr double FOOT_WIDTH = 0.04; // (m)
constexpr double FOOT_LENGTH = 0.1175; // (m)

}; // namespace DigitModified
}; // namespace RAPTOR