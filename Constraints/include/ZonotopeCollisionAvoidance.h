#ifndef ZONOTOPE_COLLISION_AVOIDANCE_H
#define ZONOTOPE_COLLISION_AVOIDANCE_H

#include "CollisionAvoidance.h"
#include <vector>

// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED
// DEPRECATED

namespace IDTO {

namespace ZonotopeParams {
    constexpr int MAX_OBSTACLE_GENERATOR_NUM = 3;
    constexpr int HYPERPLANE_NUM = MAX_OBSTACLE_GENERATOR_NUM * (MAX_OBSTACLE_GENERATOR_NUM - 1);
    constexpr int COMB_NUM = HYPERPLANE_NUM / 2;
    constexpr int VERTICES_NUM = HYPERPLANE_NUM * 4;
};

class ZonotopeCollisionAvoidance : public CollisionAvoidance {
public:
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    const int MAX_OBSTACLE_GENERATOR_NUM = ZonotopeParams::MAX_OBSTACLE_GENERATOR_NUM;
    const int HYPERPLANE_NUM = ZonotopeParams::HYPERPLANE_NUM;
    const int COMB_NUM = ZonotopeParams::COMB_NUM;
    const int VERTICES_NUM = ZonotopeParams::VERTICES_NUM;

    // Constructor
    ZonotopeCollisionAvoidance() = default;

    ZonotopeCollisionAvoidance(const Eigen::Array<Vec3, 1, Eigen::Dynamic>& zonotopeCenters_input,
                               const Eigen::Array<MatX, 1, Eigen::Dynamic>& zonotopeGenerators_input);

    // Destructor
    ~ZonotopeCollisionAvoidance() {
        delete [] A;
        delete [] b;
        delete [] v_start_idx;
        delete [] v_size;
        delete [] v;
    }

    // class methods:
    void initialize();

    void computeDistance(const Vec3& point) override;

    void computeDistance(const Vec3& point, const MatX& ppoint_pz) override;

    // class members:
        // zonotope centers and generators
    Eigen::Array<Vec3, 1, Eigen::Dynamic> zonotopeCenters;
    Eigen::Array<MatX, 1, Eigen::Dynamic> zonotopeGenerators;

        // hyperplane representation
    double* A = nullptr;
    double* b = nullptr;

        // vertices representation
    int* v_start_idx = nullptr;
    int* v_size = nullptr;
    double* v = nullptr;
};

}; // namespace IDTO

#endif // ZONOTOPE_COLLISION_AVOIDANCE_H
