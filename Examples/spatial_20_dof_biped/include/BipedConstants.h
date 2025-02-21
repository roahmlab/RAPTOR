#ifndef BIPED_CONSTANTS_H
#define BIPED_CONSTANTS_H

#include <cmath>

namespace RAPTOR {
namespace Biped {

// constants related to talos
constexpr int NUM_JOINTS = 20;
constexpr int NUM_DEPENDENT_JOINTS = 6;
constexpr int NUM_INDEPENDENT_JOINTS = 14;

const double JOINT_LIMITS_LOWER[NUM_JOINTS] = {
    -1000,           // Px
    -1000,           // Py
    -1000,           // Pz
    -1000,           // Rx
    -1000,           // Ry
    -1000,           // Rz
    -M_PI,    // stance_ab/ad
    -M_PI,     // stance hip
    -M_PI, // stance_knee (deg2rad(5))
    -M_PI,          // swing_ab/ad
    -M_PI,     // swing hip
    -M_PI, // swing_knee (deg2rad(5))
    -M_PI / 100,   // shoulder yaw
    -M_PI / 100,   // head
    -M_PI / 100,   // shoulder 1 in/out
    -M_PI / 10,    // shoulder 1 back/forth
    -M_PI / 10,    // elbow 1
    -M_PI / 100,   // shoulder 2 in/out
    -M_PI / 10,    // shoulder 2 back/forth
    -M_PI / 10     // elbow 2
};

const double JOINT_LIMITS_UPPER[NUM_JOINTS] = {
     1000,           // Px
     1000,           // Py
     1000,           // Pz
     1000,           // Rx
     1000,           // Ry
     1000,           // Rz
     M_PI,          // stance_ab/ad
     M_PI,     // stance hip
     M_PI,    // stance_knee
     M_PI,    // swing_ab/ad
     M_PI,     // swing hip
     M_PI,     // swing_knee
     M_PI / 100,   // shoulder yaw
     M_PI / 100,   // head
     M_PI / 100,   // shoulder 1 in/out
     M_PI / 10,    // shoulder 1 back/forth
    -M_PI / 30,    // elbow 1
     M_PI / 100,   // shoulder 2 in/out
     M_PI / 10,    // shoulder 2 back/forth
    -M_PI / 30     // elbow 2
};

constexpr double TORQUE_LIMITS_LOWER[NUM_INDEPENDENT_JOINTS] = {
    -300, 
    -300,
    -300,
    -300,
    -300,
    -300,
    -300,
    -300,
    -300,
    -300,
    -300,
    -300,
    -300,
    -300
};

constexpr double TORQUE_LIMITS_UPPER[NUM_INDEPENDENT_JOINTS] = {
    300,
    300,
    300,
    300,
    300,
    300,
    300,
    300,
    300,
    300,
    300,
    300,
    300,
    300
};

constexpr double MU = 0.7;
constexpr double GAMMA = 0.7;
constexpr double FOOT_WIDTH = 0.10; // (m)
constexpr double FOOT_LENGTH = 0.03; // (m)

constexpr char LEFT_FOOT_NAME[] = "left_leg_3";
constexpr char RIGHT_FOOT_NAME[] = "right_leg_3";

}; // namespace Biped
}; // namespace RAPTOR

#endif // BIPED_CONSTANTS_H