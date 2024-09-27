#include "BoxCollisionAvoidance.h"

namespace RAPTOR {

namespace Box{ 
void TensorProduct(const Eigen::Matrix3d& R, 
                   const Eigen::Array<Eigen::MatrixXd, 3, 1>& inp,
                   Eigen::Array<Eigen::MatrixXd, 3, 1>& out) {
	out(0) = R(0, 0) * inp(0) + R(0, 1) * inp(1) + R(0, 2) * inp(2);
	out(1) = R(1, 0) * inp(0) + R(1, 1) * inp(1) + R(1, 2) * inp(2);
	out(2) = R(2, 0) * inp(0) + R(2, 1) * inp(1) + R(2, 2) * inp(2);
}
}; // namespace Box

double distancePointAndLineSegment(const Eigen::Vector3d& point, 
                                   const Eigen::Vector3d& p1, 
                                   const Eigen::Vector3d& p2) {
    Eigen::Vector3d v = p2 - p1; //!< Vector representing the line segment
    Eigen::Vector3d w = point - p1; //!< Vector from p1 to the point

    double lambda = w.dot(v) / v.dot(v);
    double t = fmax(0.0, fmin(1.0, lambda));
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
	normals.resize(numObstacles, Box::HYPERPLANE_NUM);
	intercepts.resize(numObstacles, Box::HYPERPLANE_NUM);
	vertices.resize(numObstacles, Box::VERTICES_NUM);

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

Eigen::Vector3d BoxCollisionAvoidance::computeDifferenceWithCloestPoint(const Vec3& point, 
											   			  				const int obs_id,
																		double& isInside) const {
	// Compute the closest point in the local frame of the box
	Vec3 localPoint = boxR(obs_id).transpose() * (point - boxCenters[obs_id]);
	Vec3 localClosestPoint;
	isInside = -1.0;
	double distancesWithFace[3][2];

	for (int i = 0; i < 3; i++) {
		distancesWithFace[i][0] = localPoint(i) + boxSize[obs_id](i) / 2.0;
		distancesWithFace[i][1] = boxSize[obs_id](i) / 2.0 - localPoint(i);
		localClosestPoint(i) = fmax(-boxSize[obs_id](i) / 2.0, fmin(localPoint(i), boxSize[obs_id](i) / 2.0));

		if (localClosestPoint(i) != localPoint(i)) {
			isInside = 1.0;
		}
	}

	// The point is inside the box, find the closest point on the surface
	if (isInside == -1.0) {
		double minDist = 1e19;
		size_t faceIndex = 0, directionIndex = 0;
		for (int i = 0; i < 3; i++) {
			if (distancesWithFace[i][0] < minDist) {
				minDist = distancesWithFace[i][0];
				faceIndex = i;
				directionIndex = 0;
			}
			if (distancesWithFace[i][1] < minDist) {
				minDist = distancesWithFace[i][1];
				faceIndex = i;
				directionIndex = 1;
			}
		}

		localClosestPoint(faceIndex) = (directionIndex == 0) ? 
											-boxSize[obs_id](faceIndex) / 2.0 : 
											boxSize[obs_id](faceIndex) / 2.0;
	}

	// Transform the closest point from the local frame to the global frame
	Vec3 closestPoint = boxR(obs_id) * localClosestPoint + boxCenters[obs_id];

	return point - closestPoint;
}

Eigen::Vector3d BoxCollisionAvoidance::computeDifferenceWithCloestPoint(const Vec3& point, 
																		const MatX& ppoint_pz,
																		const int obs_id,
																		MatX& pdiff_pz,
																		double& isInside) const {
	pdiff_pz.resize(3, ppoint_pz.cols());	

	// Compute the closest point in the local frame of the box
	Vec3 localPoint = boxR(obs_id).transpose() * (point - boxCenters[obs_id]);
	MatX plocalPoint_pz = boxR(obs_id).transpose() * ppoint_pz;
	Vec3 localClosestPoint;
	MatX plocalClosestPoint_pz(3, ppoint_pz.cols());
	isInside = -1.0;
	double distancesWithFace[3][2];

	for (int i = 0; i < 3; i++) {
		distancesWithFace[i][0] = localPoint(i) + boxSize[obs_id](i) / 2.0;
		distancesWithFace[i][1] = boxSize[obs_id](i) / 2.0 - localPoint(i);

		if (localPoint(i) < -boxSize[obs_id](i) / 2.0) {
			localClosestPoint(i) = -boxSize[obs_id](i) / 2.0;
			isInside = 1.0;
			plocalClosestPoint_pz.row(i).setZero();
		}
		else if (localPoint(i) > boxSize[obs_id](i) / 2.0) {
			localClosestPoint(i) = boxSize[obs_id](i) / 2.0;
			isInside = 1.0;
			plocalClosestPoint_pz.row(i).setZero();
		}
		else {
			localClosestPoint(i) = localPoint(i);
			plocalClosestPoint_pz.row(i) = plocalPoint_pz.row(i);
		}
	}

	// The point is inside the box, find the closest point on the surface
	if (isInside == -1.0) {
		double minDist = 1e19;
		size_t faceIndex = 0, directionIndex = 0;
		for (int i = 0; i < 3; i++) {
			if (distancesWithFace[i][0] < minDist) {
				minDist = distancesWithFace[i][0];
				faceIndex = i;
				directionIndex = 0;
			}
			if (distancesWithFace[i][1] < minDist) {
				minDist = distancesWithFace[i][1];
				faceIndex = i;
				directionIndex = 1;
			}
		}

		localClosestPoint(faceIndex) = (directionIndex == 0) ? 
											-boxSize[obs_id](faceIndex) / 2.0 : 
											boxSize[obs_id](faceIndex) / 2.0;
		plocalClosestPoint_pz.row(faceIndex).setZero();
	}

	// Transform the closest point from the local frame to the global frame
	Vec3 closestPoint = boxR(obs_id) * localClosestPoint + boxCenters[obs_id];
	MatX pclosestPoint_pz = boxR(obs_id) * plocalClosestPoint_pz;
	
	pdiff_pz = ppoint_pz - pclosestPoint_pz;

	return point - closestPoint;											
}

Eigen::Vector3d BoxCollisionAvoidance::computeDifferenceWithCloestPoint(const Vec3& point, 
																		const MatX& ppoint_pz,
																		const Eigen::Array<MatX, 3, 1>& ppoint_pz_pz,
																		const int obs_id,
																		MatX& pdiff_pz,
																		Eigen::Array<MatX, 3, 1>& pdiff_pz_pz,
																		double& isInside) const {
	pdiff_pz.resize(3, ppoint_pz.cols());	
	for (int i = 0; i < 3; i++) {
		pdiff_pz_pz(i).resize(ppoint_pz.cols(), ppoint_pz.cols());
	}

	// Compute the closest point in the local frame of the box
	Vec3 localPoint = boxR(obs_id).transpose() * (point - boxCenters[obs_id]);
	MatX plocalPoint_pz = boxR(obs_id).transpose() * ppoint_pz;
	Eigen::Array<MatX, 3, 1> plocalPoint_pz_pz;
	Box::TensorProduct(boxR(obs_id).transpose(), ppoint_pz_pz, plocalPoint_pz_pz);
	Vec3 localClosestPoint;
	MatX plocalClosestPoint_pz(3, ppoint_pz.cols());
	Eigen::Array<MatX, 3, 1> plocalClosestPoint_pz_pz;
	isInside = -1.0;
	double distancesWithFace[3][2];

	for (int i = 0; i < 3; i++) {
		distancesWithFace[i][0] = localPoint(i) + boxSize[obs_id](i) / 2.0;
		distancesWithFace[i][1] = boxSize[obs_id](i) / 2.0 - localPoint(i);

		if (localPoint(i) < -boxSize[obs_id](i) / 2.0) {
			localClosestPoint(i) = -boxSize[obs_id](i) / 2.0;
			isInside = 1.0;
			plocalClosestPoint_pz.row(i).setZero();
			plocalClosestPoint_pz_pz(i) = MatX::Zero(ppoint_pz.cols(), ppoint_pz.cols());
		}
		else if (localPoint(i) > boxSize[obs_id](i) / 2.0) {
			localClosestPoint(i) = boxSize[obs_id](i) / 2.0;
			isInside = 1.0;
			plocalClosestPoint_pz.row(i).setZero();
			plocalClosestPoint_pz_pz(i) = MatX::Zero(ppoint_pz.cols(), ppoint_pz.cols());
		}
		else {
			localClosestPoint(i) = localPoint(i);
			plocalClosestPoint_pz.row(i) = plocalPoint_pz.row(i);
			plocalClosestPoint_pz_pz(i) = plocalPoint_pz_pz(i);
		}
	}

	// The point is inside the box, find the closest point on the surface
	if (isInside == -1.0) {
		double minDist = 1e19;
		size_t faceIndex = 0, directionIndex = 0;
		for (int i = 0; i < 3; i++) {
			if (distancesWithFace[i][0] < minDist) {
				minDist = distancesWithFace[i][0];
				faceIndex = i;
				directionIndex = 0;
			}
			if (distancesWithFace[i][1] < minDist) {
				minDist = distancesWithFace[i][1];
				faceIndex = i;
				directionIndex = 1;
			}
		}

		localClosestPoint(faceIndex) = (directionIndex == 0) ? 
											-boxSize[obs_id](faceIndex) / 2.0 : 
											boxSize[obs_id](faceIndex) / 2.0;
		plocalClosestPoint_pz.row(faceIndex).setZero();
		plocalClosestPoint_pz_pz(faceIndex).setZero();
	}

	// Transform the closest point from the local frame to the global frame
	Vec3 closestPoint = boxR(obs_id) * localClosestPoint + boxCenters[obs_id];
	MatX pclosestPoint_pz = boxR(obs_id) * plocalClosestPoint_pz;
	Eigen::Array<MatX, 3, 1> pclosestPoint_pz_pz;
	Box::TensorProduct(boxR(obs_id), plocalClosestPoint_pz_pz, pclosestPoint_pz_pz);
	
	pdiff_pz = ppoint_pz - pclosestPoint_pz;
	for (int i = 0; i < 3; i++) {
		pdiff_pz_pz(i) = ppoint_pz_pz(i) - pclosestPoint_pz_pz(i);
	}

	return point - closestPoint;
}

void BoxCollisionAvoidance::computeDistance(const Vec3& point) {
	minimumDistance = 1e19;
	Vec3 diff;
	double isInside = 1.0;

	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		diff = computeDifferenceWithCloestPoint(point, 
												obs_id, 
												isInside);

		const double diffSquared = diff.dot(diff);
		distances(obs_id) = (diffSquared > Box::SQUARE_ROOT_THRESHOLD) ?
								isInside * std::sqrt(diffSquared) : 
								0.0;

		if (distances(obs_id) < minimumDistance) {
			minimumDistance = distances(obs_id);
			minimumDistanceIndex = obs_id;
		}
	}
}

void BoxCollisionAvoidance::computeDistance(const Vec3& point, 
										    const MatX& ppoint_pz) {
	if (ppoint_pz.rows() != 3) {
		throw std::invalid_argument("ppoint_pz should have 3 rows");
	}

	minimumDistance = 1e19;
	Vec3 diff;
	MatX pdiff_pz;
	double isInside = 1.0;

	pdistances_pz.resize(numObstacles, ppoint_pz.cols());

	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		if (onlyComputeDerivativesForMinimumDistance) {
			diff = computeDifferenceWithCloestPoint(point, 
													obs_id,
													isInside);
		}
		else {
			diff = computeDifferenceWithCloestPoint(point, ppoint_pz, 
													obs_id, 
													pdiff_pz,
													isInside);
		}
		
		const double diffSquared = diff.dot(diff);
		const double distanceNorm = std::sqrt(diffSquared);
		distances(obs_id) = (diffSquared > Box::SQUARE_ROOT_THRESHOLD) ?
								isInside * distanceNorm : 
								0.0;

		if (distances(obs_id) < minimumDistance) {
			minimumDistance = distances(obs_id);
			minimumDistanceIndex = obs_id;
		}

		if (!onlyComputeDerivativesForMinimumDistance) {
			if (diffSquared > Box::SQUARE_ROOT_THRESHOLD) {
				pdistances_pz.row(obs_id) = isInside * diff.transpose() * pdiff_pz / distanceNorm;
			}
			else {
				pdistances_pz.row(obs_id).setZero();
			}
		}
	}

	if (onlyComputeDerivativesForMinimumDistance) {
		diff = computeDifferenceWithCloestPoint(point, ppoint_pz, 
												minimumDistanceIndex, 
												pdiff_pz,
												isInside);

		const double diffSquared = diff.dot(diff);
		if (diffSquared > Box::SQUARE_ROOT_THRESHOLD) {
			pdistances_pz.row(minimumDistanceIndex) = isInside * diff.transpose() * pdiff_pz / std::sqrt(diffSquared);
		}
		else {
			pdistances_pz.row(minimumDistanceIndex).setZero();
		}
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

	minimumDistance = 1e19;
	Vec3 diff;
	MatX pdiff_pz;
	Eigen::Array<MatX, 3, 1> pdiff_pz_pz;
	double isInside = 1.0;

	pdistances_pz.resize(numObstacles, ppoint_pz.cols());

	pdistances_pz_pz.resize(numObstacles);
	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		pdistances_pz_pz(obs_id).resize(ppoint_pz.cols(), ppoint_pz.cols());
	}

	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		if (onlyComputeDerivativesForMinimumDistance) {
			diff = computeDifferenceWithCloestPoint(point, 
													obs_id,
													isInside);
		}
		else {
			diff = computeDifferenceWithCloestPoint(point, ppoint_pz, ppoint_pz_pz,
													obs_id, 
													pdiff_pz, pdiff_pz_pz,
													isInside);
		}

		const double diffSquared = diff.dot(diff);
		const double distanceNorm = std::sqrt(diffSquared);
		distances(obs_id) = (diffSquared > Box::SQUARE_ROOT_THRESHOLD) ?
								isInside * distanceNorm : 
								0.0;

		if (distances(obs_id) < minimumDistance) {
			minimumDistance = distances(obs_id);
			minimumDistanceIndex = obs_id;
		}

		if (!onlyComputeDerivativesForMinimumDistance) {
			if (diffSquared > Box::SQUARE_ROOT_THRESHOLD) {
				const MatX pdiffSquare_pz = diff.transpose() * pdiff_pz;
				pdistances_pz.row(obs_id) = isInside * pdiffSquare_pz / distanceNorm;

				pdistances_pz_pz(obs_id) = pdiff_pz.transpose() * pdiff_pz / distanceNorm;
				for (int i = 0; i < 3; i++) {
					pdistances_pz_pz(obs_id) += diff(i) * pdiff_pz_pz(i) / distanceNorm;
				}
				pdistances_pz_pz(obs_id) -= pdiffSquare_pz.transpose() * pdiffSquare_pz / std::pow(distanceNorm, 3);
				pdistances_pz_pz(obs_id) *= isInside;
			}
			else {
				pdistances_pz.row(obs_id).setZero();
				pdistances_pz_pz(obs_id).setZero();
			}
		}
	}

	if (onlyComputeDerivativesForMinimumDistance) {
		diff = computeDifferenceWithCloestPoint(point, ppoint_pz, ppoint_pz_pz,
												minimumDistanceIndex, 
												pdiff_pz, pdiff_pz_pz,
												isInside);
												
		double diffSquared = diff.dot(diff);
		if (diffSquared > Box::SQUARE_ROOT_THRESHOLD) {
			double dist = std::sqrt(diffSquared);
			MatX pdiffSquare_pz = diff.transpose() * pdiff_pz;
			pdistances_pz.row(minimumDistanceIndex) = isInside * pdiffSquare_pz / dist;

			pdistances_pz_pz(minimumDistanceIndex) = pdiff_pz.transpose() * pdiff_pz / dist;
			for (int i = 0; i < 3; i++) {
				pdistances_pz_pz(minimumDistanceIndex) += diff(i) * pdiff_pz_pz(i) / dist;
			}
			pdistances_pz_pz(minimumDistanceIndex) -= pdiffSquare_pz.transpose() * pdiffSquare_pz / std::pow(diffSquared, 1.5);
			pdistances_pz_pz(minimumDistanceIndex) *= isInside;
		}
		else {
			pdistances_pz.row(minimumDistanceIndex).setZero();
			pdistances_pz_pz(minimumDistanceIndex).setZero();
		}
	}
}

}; // namespace RAPTOR
