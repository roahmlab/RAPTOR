#ifndef COLLISION_AVOIDANCE_H
#define COLLISION_AVOIDANCE_H

#include "Constraints.h"
#include "Qhull.h"

namespace IDTO {

#define MAX_OBSTACLE_NUM 40

/*
CollisionAvoidance class is one special class that does not inherit from Constraints class.
It is used to compute the distance between a point and all obstacles.
*/

class CollisionAvoidance {
public:
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    CollisionAvoidance() = default;

    // Destructor
    ~CollisionAvoidance() = default;

    // class methods:
    virtual void computeDistance(const Vec3& point) = 0;

    virtual void computeDistance(const Vec3& point, const MatX& ppoint_pz) = 0;

    // class members:
    int num_obstacles = 0;

    // compute results are stored here
    VecX g;
    MatX pg_pz;
};

}; // namespace IDTO

#endif // COLLISION_AVOIDANCE_H
