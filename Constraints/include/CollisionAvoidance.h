#ifndef COLLISION_AVOIDANCE_H
#define COLLISION_AVOIDANCE_H

#include "Constraints.h"

namespace RAPTOR {

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

    virtual void computeDistance(const Vec3& point, 
                                 const MatX& ppoint_pz) = 0;

    virtual void computeDistance(const Vec3& point, 
                                 const MatX& ppoint_pz,
                                 const Eigen::Array<MatX, 3, 1>& ppoint_pz_pz) = 0;

    // class members:
    int numObstacles = 0;

    // compute results are stored here
    VecX distances;
    MatX pdistances_pz;
    Eigen::Array<MatX, Eigen::Dynamic, 1> pdistances_pz_pz;

    double minimumDistance = 0;
    size_t minimumDistanceIndex = 0;

    bool onlyComputeDerivativesForMinimumDistance = false;
};

}; // namespace RAPTOR

#endif // COLLISION_AVOIDANCE_H
