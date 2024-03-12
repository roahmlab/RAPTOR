#ifndef BOX_COLLISION_AVOIDANCE_H
#define BOX_COLLISION_AVOIDANCE_H

#include "CollisionAvoidance.h"

namespace IDTO {

#define PROJECT_POINT_ON_FACE_THRESHOLD 1e-4

namespace Box {
    constexpr int HYPERPLANE_NUM = 6;
    constexpr int VERTICES_NUM = 8;
};

double distancePointLineSegment(const Eigen::Vector3d& point, 
                                const Eigen::Vector3d& p1, 
                                const Eigen::Vector3d& p2);

class BoxCollisionAvoidance : public CollisionAvoidance {
public:
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    const int HYPERPLANE_NUM = Box::HYPERPLANE_NUM;
    const int VERTICES_NUM = Box::VERTICES_NUM;

    // Constructor
    BoxCollisionAvoidance() = default;

    BoxCollisionAvoidance(const Eigen::Array<Vec3, 1, Eigen::Dynamic>& boxCenters_input,
                          const Eigen::Array<Vec3, 1, Eigen::Dynamic>& boxOrientation_input,
                          const Eigen::Array<Vec3, 1, Eigen::Dynamic>& boxSize_input);

    // Destructor
    ~BoxCollisionAvoidance() = default;

    // class methods:
    void initialize();

    void computeDistance(const Vec3& point) override;

    void computeDistance(const Vec3& point, const MatX& ppoint_pz) override;

    // class members:
        // box centers and generators
    Eigen::Array<Vec3, 1, Eigen::Dynamic> boxCenters;
    Eigen::Array<Vec3, 1, Eigen::Dynamic> boxOrientation;
    Eigen::Array<Vec3, 1, Eigen::Dynamic> boxSize;
    Eigen::Array<Mat3, 1, Eigen::Dynamic> boxR;

        // hyperplane representation
    Eigen::Array<Vec3, Eigen::Dynamic, Eigen::Dynamic> normals; //!< Normals of the box faces
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> intercepts; //!< Intercepts of the box faces

        // vertices representation
    Eigen::Array<Vec3, Eigen::Dynamic, Eigen::Dynamic> vertices; //!< Vertices of the box
};

}; // namespace IDTO

#endif // BOX_COLLISION_AVOIDANCE_H
