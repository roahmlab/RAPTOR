#include "BoxCollisionAvoidance.h"

namespace IDTO {

double distancePointLineSegment(const Eigen::Vector3d& point, 
                                const Eigen::Vector3d& p1, 
                                const Eigen::Vector3d& p2) {
    Eigen::Vector3d v = p2 - p1; //!< Vector representing the line segment
    Eigen::Vector3d w = point - p1; //!< Vector from p1 to the point

    double lambda = w.dot(v) / v.dot(v);
    double t = std::max(0.0, std::min(1.0, lambda));
    Eigen::Vector3d projection = p1 + t * v; //!< Projected point on the line segment

    // Compute the distance between the projected point and the given point
    return (projection - point).norm();
}

BoxCollisionAvoidance::BoxCollisionAvoidance(const Eigen::Array<Vec3, 1, Eigen::Dynamic>& boxCenters_input,
											 const Eigen::Array<Vec3, 1, Eigen::Dynamic>& boxOrientation_input,
											 const Eigen::Array<Vec3, 1, Eigen::Dynamic>& boxSize_input) :
    boxCenters(boxCenters_input), 
    boxOrientation(boxOrientation_input),
	boxSize(boxSize_input) {
	if (boxCenters.size() != boxOrientation.size() || boxCenters.size() != boxSize.size()) {
		throw std::invalid_argument("boxCenters, boxOrientation, and boxSize should have the same size");
	}

	numObstacles = boxCenters.size();
	distances.resize(numObstacles);

	initialize();
}

void BoxCollisionAvoidance::initialize() {
	// allocate memory
	boxR.resize(numObstacles);
	normals.resize(numObstacles, HYPERPLANE_NUM);
	intercepts.resize(numObstacles, HYPERPLANE_NUM);
	vertices.resize(numObstacles, VERTICES_NUM);

	for (int i = 0; i < numObstacles; i++) {
		const Vec3 half_size = boxSize[i] / 2;
		const Vec3& center = boxCenters[i];
		const Vec3& rpy = boxOrientation[i];

		Mat3 R = (Eigen::AngleAxisd(rpy[0], Vec3::UnitX())
				  * Eigen::AngleAxisd(rpy[1], Vec3::UnitY())
				  * Eigen::AngleAxisd(rpy[2], Vec3::UnitZ())).matrix();

		// initialize vertices
		vertices(i,0) = center + R * Vec3(half_size(0), half_size(1), half_size(2));
		vertices(i,1) = center + R * Vec3(half_size(0), half_size(1), -half_size(2));
		vertices(i,2) = center + R * Vec3(half_size(0), -half_size(1), half_size(2));
		vertices(i,3) = center + R * Vec3(half_size(0), -half_size(1), -half_size(2));
		vertices(i,4) = center + R * Vec3(-half_size(0), half_size(1), half_size(2));
		vertices(i,5) = center + R * Vec3(-half_size(0), half_size(1), -half_size(2));
		vertices(i,6) = center + R * Vec3(-half_size(0), -half_size(1), half_size(2));
		vertices(i,7) = center + R * Vec3(-half_size(0), -half_size(1), -half_size(2));

		// initialize hyperplanes  
		normals(i,0) = R * Vec3(1, 0, 0);
		normals(i,1) = R * Vec3(0, 1, 0);
		normals(i,2) = R * Vec3(0, 0, 1);
		normals(i,3) = R * Vec3(-1, 0, 0);
		normals(i,4) = R * Vec3(0, -1, 0);
		normals(i,5) = R * Vec3(0, 0, -1);

		intercepts(i,0) = normals(i,0).dot(center + R * Vec3(half_size(0), 0, 0));
		intercepts(i,1) = normals(i,1).dot(center + R * Vec3(0, half_size(1), 0));
		intercepts(i,2) = normals(i,2).dot(center + R * Vec3(0, 0, half_size(2)));
		intercepts(i,3) = normals(i,3).dot(center + R * Vec3(-half_size(0), 0, 0));
		intercepts(i,4) = normals(i,4).dot(center + R * Vec3(0, -half_size(1), 0));
		intercepts(i,5) = normals(i,5).dot(center + R * Vec3(0, 0, -half_size(2)));

		boxR(i) = R;
	}
}

void BoxCollisionAvoidance::computeDistance(const Vec3& point) {
	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		// Compute the closest point in the local frame of the box
		Vec3 localPoint = boxR(obs_id).transpose() * (point - boxCenters(obs_id));
		Vec3 localClosestPoint;
		for (int i = 0; i < 3; i++) {
			localClosestPoint(i) = std::max(-boxSize(obs_id)(i) / 2.0, std::min(localPoint(i), boxSize(obs_id)(i) / 2.0));
		}

		// Transform the closest point from the local frame to the global frame
		Vec3 closestPoint = boxR(obs_id) * localClosestPoint + boxCenters(obs_id);

		distances(obs_id) = (point - closestPoint).norm();
	}
}

void BoxCollisionAvoidance::computeDistance(const Vec3& point, const MatX& ppoint_pz) {
	if (ppoint_pz.rows() != 3) {
		throw std::invalid_argument("ppoint_pz should have 3 rows");
	}

	pdistances_pz.resize(numObstacles, ppoint_pz.cols());

	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		// Compute the closest point in the local frame of the box
		Vec3 localPoint = boxR(obs_id).transpose() * (point - boxCenters(obs_id));
		Vec3 localClosestPoint;
		for (int i = 0; i < 3; i++) {
			localClosestPoint(i) = std::max(-boxSize(obs_id)(i) / 2.0, std::min(localPoint(i), boxSize(obs_id)(i) / 2.0));
		}

		// Transform the closest point from the local frame to the global frame
		Vec3 closestPoint = boxR(obs_id) * localClosestPoint + boxCenters(obs_id);

		distances(obs_id) = (point - closestPoint).norm();

		if (distances(obs_id) > 1e-5) {
			pdistances_pz.row(obs_id) = (point - closestPoint).transpose() * ppoint_pz / distances(obs_id);
		}
		else {
			pdistances_pz.row(obs_id).setZero();
		}
	}
}

}; // namespace IDTO