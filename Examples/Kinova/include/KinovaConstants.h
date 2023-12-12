namespace IDTO {
namespace Kinova {

constexpr int NUM_JOINTS = 7;

constexpr double JOINT_LIMITS_LOWER[NUM_JOINTS] = {-1e19, -2.41, -1e19, -2.66, -1e19, -2.23, -1e19}; // radian

constexpr double JOINT_LIMITS_UPPER[NUM_JOINTS] = {1e19, 2.41, 1e19, 2.66, 1e19, 2.23, 1e19}; // radian

constexpr double TORQUE_LIMITS_LOWER[NUM_JOINTS] = {-56.7, -56.7, -56.7, -56.7, -29.4, -29.4, -29.4}; // N*m

constexpr double TORQUE_LIMITS_UPPER[NUM_JOINTS] = {56.7, 56.7, 56.7, 56.7, 29.4, 29.4, 29.4}; // N*m

}; // namespace Kinova
}; // namespace IDTO