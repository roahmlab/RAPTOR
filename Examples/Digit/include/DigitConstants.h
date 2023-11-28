namespace IDTO {
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
    -80,      // left-knee
    -50.3,    // left-tarsus
    -46.2755, // left-toe-A
    -180,     // left-toe-A-rod
    -10,      // left-A2
    -45.8918, // left-toe-B
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
    -44.9815, // right-toe-A
    -180,     // right-toe-A-rod
    -10,      // right-A2
    -45.5476, // right-toe-B
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
    44.9815, // left-toe-A
    180,     // left-toe-A-rod
    10,      // left-A2
    45.5476, // left-toe-B
    180,     // left-toe-B-rod
    10,      // left-B2
    34,      // left-toe-pitch
    33,      // left-toe-roll
    60,      // right-hip-roll
    40,      // right-hip-yaw
    60,      // right-hip-pitch
    180,     // right-achilles-rod
    10,      // right-ach2
    80,      // right-knee
    50.3,    // right-tarsus
    46.2755, // right-toe-A
    180,     // right-toe-A-rod
    10,      // right-A2
    45.8918, // right-toe-B
    180,     // right-toe-B-rod
    10,      // right-B2
    44,      // right-toe-pitch
    37       // right-toe-roll
};

constexpr double TORQUE_LIMITS_LOWER[NUM_INDEPENDENT_JOINTS] = {
    -126.682,   // left-hip-roll
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

constexpr double TORQUE_LIMITS_UPPER[NUM_INDEPENDENT_JOINTS] = {
    126.682,   // left-hip-roll
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

}; // namespace Digit
}; // namespace IDTO