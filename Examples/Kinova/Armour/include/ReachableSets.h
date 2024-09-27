#ifndef REACHABLE_SETS_H
#define REACHABLE_SETS_H

#include "PZDynamics.h"

namespace RAPTOR {
namespace Armour {

// Part I: Create JRS (joint trajectory reachable sets) online
void GenerateJRS(const std::shared_ptr<RobotInfo>& robotInfoPtr_,
                 std::shared_ptr<BezierCurveInterval>& trajPtr_);

// Part II: Generate link and torque PZs
void GenerateLinkAndTorquePZs(const std::shared_ptr<RobotInfo>& robotInfoPtr_,
                              std::shared_ptr<BezierCurveInterval>& trajPtr_,
                              std::shared_ptr<KinematicsDynamics>& kdPtr_);

// Part III: Compute bounds for robust input
Eigen::MatrixXd ComputeRobustInputBounds(const std::shared_ptr<RobotInfo>& robotInfoPtr_,
                                         std::shared_ptr<BezierCurveInterval>& trajPtr_,
                                         std::shared_ptr<KinematicsDynamics>& kdPtr_);

}; // namespace Armour
}; // namespace RAPTOR

#endif // REACHABLE_SETS_H