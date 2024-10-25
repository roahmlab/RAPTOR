#ifndef PZ_DYNAMICS_H
#define PZ_DYNAMICS_H

#include "ParameterizedTrajectory.h"

namespace RAPTOR {
namespace Kinova {
namespace Armour {

class KinematicsDynamics {
public:
	using Vec3 = Eigen::Vector3d;
	using VecX = Eigen::VectorXd;
	using Mat3 = Eigen::Matrix3d;
	using MatX = Eigen::MatrixXd;

	std::shared_ptr<RobotInfo> robotInfoPtr_ = nullptr;

	std::shared_ptr<BezierCurveInterval> trajPtr_ = nullptr;

	PZsparseArray com_arr;
    PZsparseArray mass_nominal_arr;
	PZsparseArray mass_uncertain_arr;
    PZsparseArray I_nominal_arr;
    PZsparseArray I_uncertain_arr;

	// sphere center PZs
	PZsparseArray sphere_centers;
	std::vector<double> sphere_radii;

	// dynamics-related PZs
    PZsparseArray torque_nom;
    PZsparseArray torque_int;

	PZsparseArray contact_force_nom;
	PZsparseArray contact_force_int;

	PZsparseArray contact_moment_nom;
	PZsparseArray contact_moment_int;

	KinematicsDynamics(const std::shared_ptr<RobotInfo>& robot_info_input,
					   const std::shared_ptr<BezierCurveInterval>& trajPtr_input);

	void reset_trajectory(const std::shared_ptr<BezierCurveInterval>& trajPtr_input);

	// generate link PZs through forward kinematics
	void fk(const size_t s_ind);

	// generate nominal torque PZs through rnea
	void rnea(const size_t s_ind,
	          const PZsparseArray& mass_arr,
			  const PZsparseArray& I_arr,
			  PZsparseArray& u,
			  PZsparseArray& contact_force,
			  PZsparseArray& contact_moment);

	void rnea_nominal(const size_t s_ind) {
		rnea(s_ind, 
			 mass_nominal_arr, I_nominal_arr, 
			 torque_nom, contact_force_nom, contact_moment_nom);
	}

	void rnea_interval(const size_t s_ind) {
		rnea(s_ind, 
			 mass_uncertain_arr, I_uncertain_arr, 
			 torque_int, contact_force_int, contact_moment_int);
	}
};

void applyOnlyOneDimension(const std::string& jointName,
						   const PZsparse& in,
						   PZsparse& out1,
						   PZsparse& out2,
						   PZsparse& out3);

void crossOnlyOneDimension(const std::string& jointName,
						   const PZsparse& v,
						   const PZsparse& in1,
						   const PZsparse& in2,
						   const PZsparse& in3,
						   PZsparse& out1,
						   PZsparse& out2,
						   PZsparse& out3);

}; // namespace Armour
}; // namespace Kinova
}; // namespace RAPTOR

#endif