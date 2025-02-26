#ifndef PZ_DYNAMICS_H
#define PZ_DYNAMICS_H

#include "pinocchio/algorithm/model.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/regressor.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/frames.hpp"

#include "ParameterizedTrajectory.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

constexpr size_t FRICTION_CONE_LINEARIZED_SIZE = 8;
constexpr size_t ZMP_LINEARIZED_SIZE = 8;
const size_t NUM_CONTACT_CONSTRAINTS = 1 + FRICTION_CONE_LINEARIZED_SIZE + ZMP_LINEARIZED_SIZE;

class PZDynamics {
public:
	using Vec3 = Eigen::Vector3d;
	using VecX = Eigen::VectorXd;
	using Mat3 = Eigen::Matrix3d;
	using MatX = Eigen::MatrixXd;

	PZDynamics() = default;

	PZDynamics(const std::shared_ptr<RobotInfo>& robotInfoPtr_input,
			   const std::shared_ptr<BezierCurveInterval>& trajPtr_input);

	~PZDynamics() = default;

	void reset_robot_info(const std::shared_ptr<RobotInfo>& robotInfoPtr_input);

	void reset_trajectory(const std::shared_ptr<BezierCurveInterval>& trajPtr_input);

	void compute();

	std::shared_ptr<RobotInfo> robotInfoPtr_ = nullptr;
	std::shared_ptr<BezierCurveInterval> trajPtr_ = nullptr;

	VecX phi;
	Eigen::Vector<Interval, Eigen::Dynamic> phi_interval;

	std::vector<pinocchio::ModelTpl<PZSparse>> model_sparses;
	std::vector<pinocchio::DataTpl<PZSparse>> data_sparses;

	std::vector<pinocchio::ModelTpl<PZSparse>> model_sparses_interval;
	std::vector<pinocchio::DataTpl<PZSparse>> data_sparses_interval;

	Eigen::Matrix<PZSparse, Eigen::Dynamic, Eigen::Dynamic> friction_PZs;
	Eigen::Matrix<PZSparse, Eigen::Dynamic, Eigen::Dynamic> zmp_PZs;

	// the radius of the torque PZs
	MatX torque_radii;

	// the radius of the sphere reachable sets
	MatX sphere_radii;

	// hyperplanes for friction cone constraints
	Eigen::Array<Eigen::Vector3d, FRICTION_CONE_LINEARIZED_SIZE, 1> S;

	// hyperplanes for ZMP constraints
	Eigen::Array<Eigen::Vector3d, ZMP_LINEARIZED_SIZE, 1> c;
	Eigen::Array<Eigen::Vector3d, ZMP_LINEARIZED_SIZE, 1> A;
};

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR

#endif