#ifndef ZONOTOPE_COLLISION_AVOIDANCE_H
#define ZONOTOPE_COLLISION_AVOIDANCE_H

#include "CollisionAvoidance.h"

namespace IDTO {

#define MAX_OBSTACLE_GENERATOR_NUM 3
#define HYPERPLANE_NUM MAX_OBSTACLE_GENERATOR_NUM * (MAX_OBSTACLE_GENERATOR_NUM - 1)
#define COMB_NUM HYPERPLANE_NUM / 2
#define VERTICES_NUM HYPERPLANE_NUM * 4

class ZonotopeCollisionAvoidance : public CollisionAvoidance {
public:
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    ZonotopeCollisionAvoidance() = default;

    ZonotopeCollisionAvoidance(const Eigen::Array<Vec3, 1, Eigen::Dynamic>& zonotopeCenters_input,
                               const Eigen::Array<MatX, 1, Eigen::Dynamic>& zonotopeGenerators_input);

    // Destructor
    ~ZonotopeCollisionAvoidance() = default;

    // class methods:
    void initialize();

    void computeDistance(const Vec3& point) override;

    void computeDistance(const Vec3& point, const MatX& ppoint_pz) override;

    // class members:
        // zonotope centers and generators
    Eigen::Array<Vec3, 1, Eigen::Dynamic> zonotopeCenters;
    Eigen::Array<MatX, 1, Eigen::Dynamic> zonotopeGenerators;

        // hyperplane representation
    double A[MAX_OBSTACLE_NUM * HYPERPLANE_NUM * 3] = {0.0};
    double b[MAX_OBSTACLE_NUM * HYPERPLANE_NUM] = {0.0};

        // vertices representation
    int v_start_idx[MAX_OBSTACLE_NUM] = {0};
    int v_size[MAX_OBSTACLE_NUM] = {0};
    double v[MAX_OBSTACLE_NUM * VERTICES_NUM * 3] = {0.0};
};

}; // namespace IDTO

#endif // ZONOTOPE_COLLISION_AVOIDANCE_H
