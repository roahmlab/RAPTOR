#include "BoxCollisionAvoidance.h"

namespace IDTO {

double distancePointAndLineSegment(const Eigen::Vector3d& point, 
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

BoxCollisionAvoidance::BoxCollisionAvoidance(const std::vector<Vec3>& boxCenters_input,
											 const std::vector<Vec3>& boxOrientation_input,
											 const std::vector<Vec3>& boxSize_input) :
    boxCenters(boxCenters_input), 
    boxOrientation(boxOrientation_input),
	boxSize(boxSize_input) {
	if (boxCenters.size() != boxOrientation.size() || 
		boxCenters.size() != boxSize.size()) {
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

	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		const Vec3 half_size = boxSize[obs_id] / 2;
		const Vec3& center = boxCenters[obs_id];
		const Vec3& rpy = boxOrientation[obs_id];

		Mat3 R = (Eigen::AngleAxisd(rpy[2], Vec3::UnitZ()) *
				  Eigen::AngleAxisd(rpy[1], Vec3::UnitY()) *
				  Eigen::AngleAxisd(rpy[0], Vec3::UnitX())).matrix();

		// initialize vertices
		vertices(obs_id, 0) = center + R * Vec3(half_size(0), half_size(1), half_size(2));
		vertices(obs_id, 1) = center + R * Vec3(half_size(0), half_size(1), -half_size(2));
		vertices(obs_id, 2) = center + R * Vec3(half_size(0), -half_size(1), half_size(2));
		vertices(obs_id, 3) = center + R * Vec3(half_size(0), -half_size(1), -half_size(2));
		vertices(obs_id, 4) = center + R * Vec3(-half_size(0), half_size(1), half_size(2));
		vertices(obs_id, 5) = center + R * Vec3(-half_size(0), half_size(1), -half_size(2));
		vertices(obs_id, 6) = center + R * Vec3(-half_size(0), -half_size(1), half_size(2));
		vertices(obs_id, 7) = center + R * Vec3(-half_size(0), -half_size(1), -half_size(2));

		// initialize hyperplanes  
		normals(obs_id, 0) = R * Vec3(1, 0, 0);
		normals(obs_id, 1) = R * Vec3(0, 1, 0);
		normals(obs_id, 2) = R * Vec3(0, 0, 1);
		normals(obs_id, 3) = R * Vec3(-1, 0, 0);
		normals(obs_id, 4) = R * Vec3(0, -1, 0);
		normals(obs_id, 5) = R * Vec3(0, 0, -1);

		intercepts(obs_id, 0) = normals(obs_id, 0).dot(center + R * Vec3(half_size(0), 0, 0));
		intercepts(obs_id, 1) = normals(obs_id, 1).dot(center + R * Vec3(0, half_size(1), 0));
		intercepts(obs_id, 2) = normals(obs_id, 2).dot(center + R * Vec3(0, 0, half_size(2)));
		intercepts(obs_id, 3) = normals(obs_id, 3).dot(center + R * Vec3(-half_size(0), 0, 0));
		intercepts(obs_id, 4) = normals(obs_id, 4).dot(center + R * Vec3(0, -half_size(1), 0));
		intercepts(obs_id, 5) = normals(obs_id, 5).dot(center + R * Vec3(0, 0, -half_size(2)));

		boxR(obs_id) = R;
	}
}

Eigen::Vector3d BoxCollisionAvoidance::computeCloestPoint(const Vec3& point, 
                            							  const int obs_id) {
	// Compute the closest point in the local frame of the box
	Vec3 localPoint = boxR(obs_id).transpose() * (point - boxCenters[obs_id]);
	Vec3 localClosestPoint;
	for (int i = 0; i < 3; i++) {
		localClosestPoint(i) = std::max(-boxSize[obs_id](i) / 2.0, std::min(localPoint(i), boxSize[obs_id](i) / 2.0));
	}

	// Transform the closest point from the local frame to the global frame
	Vec3 closestPoint = boxR(obs_id) * localClosestPoint + boxCenters[obs_id];		

	return closestPoint;					
}

void BoxCollisionAvoidance::computeDistance(const Vec3& point) {
	minimumDistance = std::numeric_limits<double>::max();

	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		Vec3 closestPoint = computeCloestPoint(point, obs_id);

		Vec3 diff = point - closestPoint;
		// distances(obs_id) = diff.norm();
		distances(obs_id) = 0.5 * diff.dot(diff);

		if (distances(obs_id) < minimumDistance) {
			minimumDistance = distances(obs_id);
			minimumDistanceIndex = obs_id;
			minimumDistanceCloestPoint = closestPoint;
		}
	}
}

void BoxCollisionAvoidance::computeDistance(const Vec3& point, 
										    const MatX& ppoint_pz) {
	if (ppoint_pz.rows() != 3) {
		throw std::invalid_argument("ppoint_pz should have 3 rows");
	}

	minimumDistance = std::numeric_limits<double>::max();

	pdistances_pz.resize(numObstacles, ppoint_pz.cols());

	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		Vec3 closestPoint = computeCloestPoint(point, obs_id);
		Vec3 diff = point - closestPoint;

		// distances(obs_id) = diff.norm();
		distances(obs_id) = diff.dot(diff);

		if (distances(obs_id) < minimumDistance) {
			minimumDistance = distances(obs_id);
			minimumDistanceIndex = obs_id;
		}

		if (!onlyComputeDerivativesForMinimumDistance) {
			// if (distances(obs_id) > 1e-5) {
			// 	pdistances_pz.row(obs_id) = (point - closestPoint).transpose() * ppoint_pz / distances(obs_id);
			// }
			// else {
			// 	pdistances_pz.row(obs_id).setZero();
			// }
			pdistances_pz.row(obs_id) = diff.transpose() * ppoint_pz;
		}
	}

	if (onlyComputeDerivativesForMinimumDistance) {
		Vec3 closestPoint = computeCloestPoint(point, minimumDistanceIndex);
		Vec3 diff = point - closestPoint;

		// if (minimumDistance > 1e-5) {
		// 	pdistances_pz.row(minimumDistanceIndex) = (point - closestPoint).transpose() * ppoint_pz / minimumDistance;
		// }
		// else {
		// 	pdistances_pz.row(minimumDistanceIndex).setZero();
		// }	
		pdistances_pz.row(minimumDistanceIndex) = diff.transpose() * ppoint_pz;
	}
}

void BoxCollisionAvoidance::computeDistance(const Vec3& point, 
											const MatX& ppoint_pz,
											const Eigen::Array<MatX, 3, 1>& ppoint_pz_pz) {
	if (ppoint_pz.rows() != 3) {
		throw std::invalid_argument("ppoint_pz should have 3 rows");
	}

    if (ppoint_pz_pz(0).rows() != ppoint_pz.cols() || 
		ppoint_pz_pz(0).cols() != ppoint_pz.cols()) {
		throw std::invalid_argument("ppoint_pz_pz should have the same number of columns as ppoint_pz");
	}
	if (ppoint_pz_pz(1).rows() != ppoint_pz.cols() || 
		ppoint_pz_pz(1).cols() != ppoint_pz.cols()) {
		throw std::invalid_argument("ppoint_pz_pz should have the same number of columns as ppoint_pz");
	}
	if (ppoint_pz_pz(2).rows() != ppoint_pz.cols() || 
		ppoint_pz_pz(2).cols() != ppoint_pz.cols()) {
		throw std::invalid_argument("ppoint_pz_pz should have the same number of columns as ppoint_pz");
	}

	minimumDistance = std::numeric_limits<double>::max();

	pdistances_pz.resize(numObstacles, ppoint_pz.cols());

	pdistances_pz_pz.resize(numObstacles);
	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		pdistances_pz_pz(obs_id).resize(ppoint_pz.cols(), ppoint_pz.cols());
	}

	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		Vec3 closestPoint = computeCloestPoint(point, obs_id);

		Vec3 diff = point - closestPoint;
		// distances(obs_id) = diff.norm();
		distances(obs_id) = diff.dot(diff);

		if (distances(obs_id) < minimumDistance) {
			minimumDistance = distances(obs_id);
			minimumDistanceIndex = obs_id;
		}

		if (!onlyComputeDerivativesForMinimumDistance) {
			// if (distances(obs_id) > 1e-5) {
			// 	pdistances_pz.row(obs_id) = (point - closestPoint).transpose() * ppoint_pz / distances(obs_id);
			// }
			// else {
			// 	pdistances_pz.row(obs_id).setZero();
			// }
			pdistances_pz.row(obs_id) = diff.transpose() * ppoint_pz;

			pdistances_pz_pz(obs_id) = ppoint_pz.transpose() * ppoint_pz;
			for (int i = 0; i < 3; i++) {
				pdistances_pz_pz(obs_id) += diff(i) * ppoint_pz_pz(i);
			}
		}
	}

	if (onlyComputeDerivativesForMinimumDistance) {
		Vec3 closestPoint = computeCloestPoint(point, minimumDistanceIndex);
		Vec3 diff = point - closestPoint;

		// if (minimumDistance > 1e-5) {
		// 	pdistances_pz.row(minimumDistanceIndex) = (point - closestPoint).transpose() * ppoint_pz / minimumDistance;
		// }
		// else {
		// 	pdistances_pz.row(minimumDistanceIndex).setZero();
		// }	
		pdistances_pz.row(minimumDistanceIndex) = diff.transpose() * ppoint_pz;

		pdistances_pz_pz(minimumDistanceIndex) = ppoint_pz.transpose() * ppoint_pz;
		for (int i = 0; i < 3; i++) {
			pdistances_pz_pz(minimumDistanceIndex) += diff(i) * ppoint_pz_pz(i);
		}
	}
}

}; // namespace IDTO