#include "ZonotopeCollisionAvoidance.h"

namespace IDTO {

ZonotopeCollisionAvoidance::ZonotopeCollisionAvoidance(const Eigen::Array<Vec3, 1, Eigen::Dynamic>& zonotopeCenters_input,
                                                       const Eigen::Array<MatX, 1, Eigen::Dynamic>& zonotopeGenerators_input) :
    zonotopeCenters(zonotopeCenters_input), 
    zonotopeGenerators(zonotopeGenerators_input) {
	for (int i = 0; i < zonotopeGenerators.size(); i++) {
		auto it = zonotopeGenerators[i];
		if (it.rows() != 3 || it.cols() != MAX_OBSTACLE_GENERATOR_NUM) {
			throw std::invalid_argument("Zonotope generator matrix should be 3 x " + std::to_string(MAX_OBSTACLE_GENERATOR_NUM));
		}
	}

	numObstacles = zonotopeCenters.size();
	distances.resize(numObstacles);

	initialize();
}

void ZonotopeCollisionAvoidance::initialize() {
    // combination index
	int combA[COMB_NUM], combB[COMB_NUM];
	int a_id = 0, b_id = 1;
	for (int i = 0; i < COMB_NUM; i++) {
		combA[i] = a_id;
		combB[i] = b_id;

		if (b_id < MAX_OBSTACLE_GENERATOR_NUM - 1) {
			b_id++;
		}
		else {
			a_id++;
			b_id = a_id + 1;
		}
	}

	// compute hyperplane representation
	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		for (int plane_id = 0; plane_id < COMB_NUM; plane_id++) {
			int a_id = combA[plane_id];
			int b_id = combB[plane_id];

            const Vec3& c = zonotopeCenters[obs_id];
            const MatX& G = zonotopeGenerators[obs_id];

            double generator_cross[3] = {G(1, a_id) * G(2, b_id) - G(2, a_id) * G(1, b_id),
                                         G(2, a_id) * G(0, b_id) - G(0, a_id) * G(2, b_id),
                                         G(0, a_id) * G(1, b_id) - G(1, a_id) * G(0, b_id)};

			double norm_cross = sqrt(pow(generator_cross[0], 2) + 
									 pow(generator_cross[1], 2) + 
									 pow(generator_cross[2], 2));

			double C[3] = {0.0};

			if (norm_cross > 0) {
				C[0] = generator_cross[0] / norm_cross;
				C[1] = generator_cross[1] / norm_cross;
				C[2] = generator_cross[2] / norm_cross;
			}

			// A = C
			A[(obs_id * HYPERPLANE_NUM + plane_id) * 3 + 0] = C[0];
			A[(obs_id * HYPERPLANE_NUM + plane_id) * 3 + 1] = C[1];
			A[(obs_id * HYPERPLANE_NUM + plane_id) * 3 + 2] = C[2];
			A[(obs_id * HYPERPLANE_NUM + plane_id + COMB_NUM) * 3 + 0] = -C[0];
			A[(obs_id * HYPERPLANE_NUM + plane_id + COMB_NUM) * 3 + 1] = -C[1];
			A[(obs_id * HYPERPLANE_NUM + plane_id + COMB_NUM) * 3 + 2] = -C[2];

			// d = C * c
			double d = C[0] * c[0] + C[1] * c[1] + C[2] * c[2];

			// delta = sum(abs(C * G))
			double delta = 0.0;
			for (int j = 0; j < MAX_OBSTACLE_GENERATOR_NUM; j++) {
				delta += fabs(C[0] * G(0, j) + C[1] * G(1, j) + C[2] * G(2, j));
			}

			b[(obs_id * HYPERPLANE_NUM + plane_id)] = d + delta;
			b[(obs_id * HYPERPLANE_NUM + plane_id + COMB_NUM)] = -d + delta;
		}
	}

	// compute vertices
	int dim = 3;
	int v_cursor = 0;
	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		const Vec3& c = zonotopeCenters[obs_id];
		const MatX& G = zonotopeGenerators[obs_id];
		std::vector<double> vertices;
		vertices.push_back(c[0]);
		vertices.push_back(c[1]);
		vertices.push_back(c[2]);
		for (int i = 0; i < MAX_OBSTACLE_GENERATOR_NUM; ++i) {
			std::vector<double> vertices_new;
			for (int j = 0; j < vertices.size() / dim; ++j) {
				for (int k = 0; k < dim; ++k) {
					vertices_new.push_back(vertices[j * dim + k] + G(k, i));
				}
				for (int k = 0; k < dim; ++k) {
					vertices_new.push_back(vertices[j * dim + k] - G(k, i));
				}
			}
			vertices = vertices_new;
			if (i >= dim) {
				vertices_new.clear();
				try {
					orgQhull::Qhull qhull("", dim, vertices.size() / 3, vertices.data(), "");
					for(const orgQhull::QhullVertex &v: qhull.vertexList()) {
						const orgQhull::QhullPoint &qhullPt = v.point();
						auto coords = qhullPt.coordinates();
						vertices_new.push_back(coords[0]);
						vertices_new.push_back(coords[1]);
						vertices_new.push_back(coords[2]);
					}
				} 
				catch (...) {
					std::cerr << "Qhull failed !" << std::endl;
					vertices_new = vertices;
				}
				vertices = vertices_new;
			}
		}
		if (obs_id == 0) {
			v_start_idx[obs_id] = 0;
		}
		else {
			v_start_idx[obs_id] = v_start_idx[obs_id - 1] + v_size[obs_id - 1];
		}
		
		v_size[obs_id] = vertices.size() / 3;

		for (auto it : vertices) {
			v[v_cursor++] = it;
		}
	}
}

void ZonotopeCollisionAvoidance::computeDistance(const Vec3& point) {
	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		// compute distances to hyperplanes
		double Apprime_minus_b[HYPERPLANE_NUM];
		double perpendicular_distances = 1e19;
		int perpendicular_distances_id = 0;
		double projected_points[3] = {0.0};

		for (int p_id = 0; p_id < HYPERPLANE_NUM; p_id++) {
			const double A_elt_x = A[(obs_id * HYPERPLANE_NUM + p_id) * 3 + 0];
			const double A_elt_y = A[(obs_id * HYPERPLANE_NUM + p_id) * 3 + 1];
			const double A_elt_z = A[(obs_id * HYPERPLANE_NUM + p_id) * 3 + 2];
			const double b_elt = b[obs_id * HYPERPLANE_NUM + p_id];

			const double A_elt_norm = fabs(A_elt_x) + fabs(A_elt_y) + fabs(A_elt_z);

			if (A_elt_norm > 0) {
				Apprime_minus_b[p_id] = A_elt_x * point(0) + 
										A_elt_y * point(1) + 
										A_elt_z * point(2) - b_elt;
			}
			else {
				Apprime_minus_b[p_id] = -1e19;
			}
		}

		// find the maximum distance to hyperplanes
		double max_elt = -1e19;
		int max_id = 0;

		for (int i = 0; i < HYPERPLANE_NUM; i++) {
			if (Apprime_minus_b[i] > max_elt) {
				max_elt = Apprime_minus_b[i];
				max_id = i;
			}
		}
		
		if (max_elt <= 0) {
			// if the distance is negative, then the sphere is already in collision and directly assign the distance
			distances(obs_id) = max_elt;
			break;
		}
		else {
			// if the distance is positive, need to compute the distance to the face
			perpendicular_distances = max_elt;
			perpendicular_distances_id = max_id;

			// compute the projected point
			const double* max_A_elt = A + (obs_id * HYPERPLANE_NUM + max_id) * 3;
			projected_points[0] = point(0) - max_elt * max_A_elt[0];
			projected_points[1] = point(1) - max_elt * max_A_elt[1];
			projected_points[2] = point(2) - max_elt * max_A_elt[2];
		}


		// determine whether the projected point is on the face
		for (int p_id = 0; p_id < HYPERPLANE_NUM; p_id++) {
			const double A_elt_x = A[(obs_id * HYPERPLANE_NUM + p_id) * 3 + 0];
			const double A_elt_y = A[(obs_id * HYPERPLANE_NUM + p_id) * 3 + 1];
			const double A_elt_z = A[(obs_id * HYPERPLANE_NUM + p_id) * 3 + 2];
			const double b_elt = b[obs_id * HYPERPLANE_NUM + p_id];

			const double A_elt_norm = fabs(A_elt_x) + fabs(A_elt_y) + fabs(A_elt_z);
			if (A_elt_norm > 0) {
				Apprime_minus_b[p_id] = A_elt_x * projected_points[0] + 
										A_elt_y * projected_points[1] + 
										A_elt_z * projected_points[2] - b_elt;
			}
			else {
				Apprime_minus_b[p_id] = -1e19;
			}
		}

		// check if the projected point is on the face
		bool on_zonotope = true;
		for (int i = 0; i < HYPERPLANE_NUM; i++) {
			if (Apprime_minus_b[i] > 1e-4) {
				on_zonotope = false;
				break;
			}
		}
		
		if (on_zonotope) {
			// if projected point is on the face, then the distance is the perpendicular distance
			const double max_elt = perpendicular_distances;
			const int max_id = perpendicular_distances_id;

			distances(obs_id) = max_elt;
		}
		else {
			// if the projected point is not on the face, then the distance is the minimum distance to the vertices corresponding to the face
			const int selected_v_start_idx = v_start_idx[obs_id];
			const int selected_v_size = v_size[obs_id];
			double min_distance = 1e19;
			int min_v_id = 0;

			for (int i = 0; i < selected_v_size; i++) {
				const double v_x = v[(selected_v_start_idx + i) * 3 + 0];
				const double v_y = v[(selected_v_start_idx + i) * 3 + 1];
				const double v_z = v[(selected_v_start_idx + i) * 3 + 2];

				const double distance = sqrt(pow(point(0) - v_x, 2) + 
											 pow(point(1) - v_y, 2) + 
											 pow(point(2) - v_z, 2));

				if (distance < min_distance) {
					min_distance = distance;
					min_v_id = i;
				}
			}

			distances(obs_id) = min_distance;
		}
	}
}

void ZonotopeCollisionAvoidance::computeDistance(const Vec3& point, const MatX& ppoint_pz) {
	if (ppoint_pz.rows() != 3) {
		throw std::invalid_argument("ppoint_pz should have 3 rows");
	}

	pdistances_pz.resize(numObstacles, ppoint_pz.cols());

	for (int obs_id = 0; obs_id < numObstacles; obs_id++) {
		// compute distances to hyperplanes
		double Apprime_minus_b[HYPERPLANE_NUM];
		double perpendicular_distances = 1e19;
		int perpendicular_distances_id = 0;
		double projected_points[3] = {0.0};

		for (int p_id = 0; p_id < HYPERPLANE_NUM; p_id++) {
			const double A_elt_x = A[(obs_id * HYPERPLANE_NUM + p_id) * 3 + 0];
			const double A_elt_y = A[(obs_id * HYPERPLANE_NUM + p_id) * 3 + 1];
			const double A_elt_z = A[(obs_id * HYPERPLANE_NUM + p_id) * 3 + 2];
			const double b_elt = b[obs_id * HYPERPLANE_NUM + p_id];

			const double A_elt_norm = fabs(A_elt_x) + fabs(A_elt_y) + fabs(A_elt_z);

			if (A_elt_norm > 0) {
				Apprime_minus_b[p_id] = A_elt_x * point(0) + 
										A_elt_y * point(1) + 
										A_elt_z * point(2) - b_elt;
			}
			else {
				Apprime_minus_b[p_id] = -1e19;
			}
		}

		// find the maximum distance to hyperplanes
		double max_elt = -1e19;
		int max_id = 0;

		for (int i = 0; i < HYPERPLANE_NUM; i++) {
			if (Apprime_minus_b[i] > max_elt) {
				max_elt = Apprime_minus_b[i];
				max_id = i;
			}
		}
		
		if (max_elt <= 0) {
			// if the distance is negative, then the sphere is already in collision and directly assign the distance
			distances(obs_id) = max_elt;

			const double* max_A_elt = A + (obs_id * HYPERPLANE_NUM + max_id) * 3;

			pdistances_pz.row(obs_id) = max_A_elt[0] * ppoint_pz.row(0) +
									    max_A_elt[1] * ppoint_pz.row(1) +
										max_A_elt[2] * ppoint_pz.row(2);

			break;
		}
		else {
			// if the distance is positive, need to compute the distance to the face
			perpendicular_distances = max_elt;
			perpendicular_distances_id = max_id;

			// compute the projected point
			const double* max_A_elt = A + (obs_id * HYPERPLANE_NUM + max_id) * 3;
			projected_points[0] = point(0) - max_elt * max_A_elt[0];
			projected_points[1] = point(1) - max_elt * max_A_elt[1];
			projected_points[2] = point(2) - max_elt * max_A_elt[2];
		}


		// determine whether the projected point is on the face
		for (int p_id = 0; p_id < HYPERPLANE_NUM; p_id++) {
			const double A_elt_x = A[(obs_id * HYPERPLANE_NUM + p_id) * 3 + 0];
			const double A_elt_y = A[(obs_id * HYPERPLANE_NUM + p_id) * 3 + 1];
			const double A_elt_z = A[(obs_id * HYPERPLANE_NUM + p_id) * 3 + 2];
			const double b_elt = b[obs_id * HYPERPLANE_NUM + p_id];

			const double A_elt_norm = fabs(A_elt_x) + fabs(A_elt_y) + fabs(A_elt_z);
			if (A_elt_norm > 0) {
				Apprime_minus_b[p_id] = A_elt_x * projected_points[0] + 
										A_elt_y * projected_points[1] + 
										A_elt_z * projected_points[2] - b_elt;
			}
			else {
				Apprime_minus_b[p_id] = -1e19;
			}
		}

		// check if the projected point is on the face
		bool on_zonotope = true;
		for (int i = 0; i < HYPERPLANE_NUM; i++) {
			if (Apprime_minus_b[i] > 1e-4) {
				on_zonotope = false;
				break;
			}
		}
		
		if (on_zonotope) {
			// if projected point is on the face, then the distance is the perpendicular distance
			const double max_elt = perpendicular_distances;
			const int max_id = perpendicular_distances_id;

			distances(obs_id) = max_elt;

			const double* max_A_elt = A + (obs_id * HYPERPLANE_NUM + max_id) * 3;

			pdistances_pz.row(obs_id) = max_A_elt[0] * ppoint_pz.row(0) +
									    max_A_elt[1] * ppoint_pz.row(1) +
										max_A_elt[2] * ppoint_pz.row(2);
		}
		else {
			// if the projected point is not on the face, then the distance is the minimum distance to the vertices corresponding to the face
			const int selected_v_start_idx = v_start_idx[obs_id];
			const int selected_v_size = v_size[obs_id];
			double min_distance = 1e19;
			int min_v_id = 0;

			for (int i = 0; i < selected_v_size; i++) {
				const double v_x = v[(selected_v_start_idx + i) * 3 + 0];
				const double v_y = v[(selected_v_start_idx + i) * 3 + 1];
				const double v_z = v[(selected_v_start_idx + i) * 3 + 2];

				const double distance = sqrt(pow(point(0) - v_x, 2) + 
											 pow(point(1) - v_y, 2) + 
											 pow(point(2) - v_z, 2));

				if (distance < min_distance) {
					min_distance = distance;
					min_v_id = i;
				}
			}

			distances(obs_id) = min_distance;

			const double v_x = v[(selected_v_start_idx + min_v_id) * 3 + 0];
			const double v_y = v[(selected_v_start_idx + min_v_id) * 3 + 1];
			const double v_z = v[(selected_v_start_idx + min_v_id) * 3 + 2];

			pdistances_pz.row(obs_id) = ((point(0) - v_x) * ppoint_pz.row(0) +
							 			 (point(1) - v_y) * ppoint_pz.row(1) +
							 			 (point(2) - v_z) * ppoint_pz.row(2)) / min_distance;
		}
	}
}

}; // namespace IDTO