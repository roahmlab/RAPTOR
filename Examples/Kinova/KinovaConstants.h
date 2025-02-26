namespace RAPTOR {
namespace Kinova {

// You can find the following information in Kinova gen3 official documentation
// https://artifactory.kinovaapps.com/ui/api/v1/download?repoKey=generic-documentation-public&path=Documentation%252FGen3%252FTechnical%2520documentation%252FUser%2520Guide%252FEN-UG-014-Gen3-Ultra-lightweight-user-guide-r9.1.pdf

constexpr int NUM_JOINTS = 7;

 // [-1e19, 1e19] means it is a continuous joint (360 degree)
 // so we don't need to set the joint limits
constexpr double JOINT_LIMITS_LOWER[NUM_JOINTS] = {-1e19, -128.9, -1e19, -147.8, -1e19, -120.3, -1e19}; // degree
                                               
constexpr double JOINT_LIMITS_UPPER[NUM_JOINTS] = {1e19, 128.9, 1e19, 147.8, 1e19, 120.3, 1e19}; // degree

constexpr double VELOCITY_LIMITS_LOWER[NUM_JOINTS] = {-100, -100, -100, -100, -180, -180, -180}; // degree/s

constexpr double VELOCITY_LIMITS_UPPER[NUM_JOINTS] = {100, 100, 100, 100, 180, 180, 180}; // degree/s

constexpr double ACCELERATION_LIMITS_LOWER[NUM_JOINTS] = {-297.94, -297.94, -297.94, -297.94, -572.95, -572.95, -572.95}; // degree/s^2

constexpr double ACCELERATION_LIMITS_UPPER[NUM_JOINTS] = {297.94, 297.94, 297.94, 297.94, 572.95, 572.95, 572.95}; // degree/s^2

// constexpr double TORQUE_LIMITS_LOWER[NUM_JOINTS] = {-56.7, -56.7, -56.7, -56.7, -29.4, -29.4, -29.4}; // N*m, warning
constexpr double TORQUE_LIMITS_LOWER[NUM_JOINTS] = {-90.0, -90.0, -90.0, -90.0, -45.0, -45.0, -45.0}; // N*m, 90% of hard limits

// constexpr double TORQUE_LIMITS_UPPER[NUM_JOINTS] = {56.7, 56.7, 56.7, 56.7, 29.4, 29.4, 29.4}; // N*m, warning
constexpr double TORQUE_LIMITS_UPPER[NUM_JOINTS] = {90.0, 90.0, 90.0, 90.0, 45.0, 45.0, 45.0}; // N*m, 90% of hard limits

constexpr double GRAVITY = -9.81; // m/s^2

}; // namespace Kinova
}; // namespace RAPTOR