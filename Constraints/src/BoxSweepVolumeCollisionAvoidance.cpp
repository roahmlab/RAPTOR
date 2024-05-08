#include "BoxSweepVolumeCollisionAvoidance.h"

namespace IDTO {

BoxSweepVolumeCollisionAvoidance::BoxSweepVolumeCollisionAvoidance(const Vec3& targetBoxCenter,
                                                                   const Vec3& targetBoxOrientation,
                                                                   const Vec3& targetBoxSize,
                                                                   const Eigen::Array<Vec3, 1, Eigen::Dynamic>& boxCenters_input,
                                                                   const Eigen::Array<Vec3, 1, Eigen::Dynamic>& boxOrientation_input,
                                                                   const Eigen::Array<Vec3, 1, Eigen::Dynamic>& boxSize_input) :
    targetBoxCenter(targetBoxCenter),
    targetBoxOrientation(targetBoxOrientation),
    targetBoxSize(targetBoxSize),
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

void BoxSweepVolumeCollisionAvoidance::initialize() {
	// allocate memory
	boxR.resize(numObstacles);
	normals.resize(numObstacles, HYPERPLANE_NUM);
	intercepts.resize(numObstacles, HYPERPLANE_NUM);
    nominators.resize(HYPERPLANE_NUM);
    denominators.resize(HYPERPLANE_NUM);
    lambda_lb.resize(HYPERPLANE_NUM);
    lambda_ub.resize(HYPERPLANE_NUM);

    Mat3 R;
    R = (Eigen::AngleAxisd(targetBoxOrientation[0], Vec3::UnitX())
         * Eigen::AngleAxisd(targetBoxOrientation[1], Vec3::UnitY())
         * Eigen::AngleAxisd(targetBoxOrientation[2], Vec3::UnitZ())).matrix();

    // the generators of two boxes combined
    Eigen::Array<Vec3, 1, Eigen::Dynamic> all_generators(GENERATOR_NUM);

    // initialize generators for target box
    all_generators(0) = R * Vec3(targetBoxSize[0] / 2, 0, 0);
    all_generators(1) = R * Vec3(0, targetBoxSize[1] / 2, 0);
    all_generators(2) = R * Vec3(0, 0, targetBoxSize[2] / 2);

	for (size_t i = 0; i < numObstacles; i++) {
		const Vec3 half_size = boxSize[i] / 2;
		const Vec3& center = boxCenters[i];
		const Vec3& rpy = boxOrientation[i];

		R = (Eigen::AngleAxisd(rpy[0], Vec3::UnitX())
             * Eigen::AngleAxisd(rpy[1], Vec3::UnitY())
             * Eigen::AngleAxisd(rpy[2], Vec3::UnitZ())).matrix();
        
        // initialize generators for the box
        all_generators(3) = R * Vec3(half_size(0), 0, 0);
        all_generators(4) = R * Vec3(0, half_size(1), 0);
        all_generators(5) = R * Vec3(0, 0, half_size(2));

        // go through a combination of any two generators in all_generators
        // to initialize the hyperplanes of the buffered obstacle
        size_t index = 0;
        for (size_t j = 0; j < all_generators.size(); j++) {
            for (size_t k = j + 1; k < all_generators.size(); k++) {
                Vec3 normal = all_generators(j).cross(all_generators(k));
                if (normal.norm() < 1e-6) {
                    normal.setZero();
                }
                else {
                    normal.normalize();
                }

                double d = normal.dot(center);
                double delta = 0;

                for (size_t l = 0; l < all_generators.size(); l++) {
                    delta += fabs(normal.dot(all_generators(l)));
                }

                normals(i, 2 * index) = normal;
                intercepts(i, 2 * index) = d + delta;

                normals(i, 2 * index + 1) = -normal;
                intercepts(i, 2 * index + 1) = -d + delta;

                index++;
            }
        }

		boxR(i) = R;
	}
}

void BoxSweepVolumeCollisionAvoidance::computeDistance(const Vec3& point) {
    for (size_t o = 0; o < numObstacles; o++) {
        for (size_t i = 0; i < HYPERPLANE_NUM; i++) {
            denominators(i) = normals(o, i).dot(point - targetBoxCenter);
            nominators(i) = intercepts(o, i) - normals(o, i).dot(targetBoxCenter);

            if (fabs(denominators(i)) < 1e-6) { // straight line segment is in parallel with this pair of hyperplanes 
                if (nominators(i) < 0) { // already impossible to collide
                    lambda_lb(i) = 1e19;
                    lambda_ub(i) = -1e19;
                }
                else { // still possible to collide
                    lambda_lb(i) = -1e19;
                    lambda_ub(i) = 1e19;
                }
            }
            else {
                if (denominators(i) > 0) {
                    lambda_lb(i) = -1e19;
                    lambda_ub(i) = nominators(i) / denominators(i);
                }
                else {
                    lambda_lb(i) = nominators(i) / denominators(i);
                    lambda_ub(i) = 1e19;
                }
            }
        }

        Eigen::Index lambda_lb_max_index = 0;
        double lambda_lb_max = lambda_lb.maxCoeff(&lambda_lb_max_index);
        Eigen::Index lambda_ub_min_index = 0;
        double lambda_ub_min = lambda_ub.minCoeff(&lambda_ub_min_index);


    }
}

void BoxSweepVolumeCollisionAvoidance::computeDistance(const Vec3& point, const MatX& ppoint_pz) {

}

}; // namespace IDTO