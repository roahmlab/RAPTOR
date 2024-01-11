#include "ZonotopeCollisionAvoidance.h"

namespace IDTO {

ZonotopeCollisionAvoidance::ZonotopeCollisionAvoidance(const Eigen::Array<Vec3, 1, Eigen::Dynamic>& zonotopeCenters_input,
                                                       const Eigen::Array<MatX, 1, Eigen::Dynamic>& zonotopeGenerators_input) :
    zonotopeCenters(zonotopeCenters_input), 
    zonotopeGenerators(zonotopeGenerators_input) {
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
	for (int obs_id = 0; obs_id < num_obstacles; obs_id++) {
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
	for (int obs_id = 0; obs_id < num_obstacles; obs_id++) {
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
						auto coords = qhullPt.coordinates(); // list of doubles
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
	// TODO
}

void ZonotopeCollisionAvoidance::computeDistance(const Vec3& point, const MatX& ppoint_pz) {
	// TODO
}

}; // namespace IDTO