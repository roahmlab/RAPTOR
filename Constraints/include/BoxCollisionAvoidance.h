#ifndef BOX_COLLISION_AVOIDANCE_H
#define BOX_COLLISION_AVOIDANCE_H

#include "CollisionAvoidance.h"
#include <vector>

namespace RAPTOR {

#define PROJECT_POINT_ON_FACE_THRESHOLD 1e-4

namespace Box {
constexpr double SQUARE_ROOT_THRESHOLD = 1e-12;

constexpr int HYPERPLANE_NUM = 6; // 3 * (3 - 1)
constexpr int VERTICES_NUM = 8;

void TensorProduct(const Eigen::Matrix3d& R, 
                   const Eigen::Array<Eigen::MatrixXd, 3, 1>& inp,
                   Eigen::Array<Eigen::MatrixXd, 3, 1>& out);
};

double distancePointAndLineSegment(const Eigen::Vector3d& point, 
                                  const Eigen::Vector3d& p1, 
                                  const Eigen::Vector3d& p2);

class BoxCollisionAvoidance : public CollisionAvoidance {
public:
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    // Constructor
    BoxCollisionAvoidance() = default;

    BoxCollisionAvoidance(const std::vector<Vec3>& boxCenters_input,
                          const std::vector<Vec3>& boxOrientation_input,
                          const std::vector<Vec3>& boxSize_input);

    // Destructor
    ~BoxCollisionAvoidance() = default;

    // class methods:
    void initialize();

    Vec3 computeDifferenceWithCloestPoint(const Vec3& point, 
                                          const int obs_id,
                                          double& isInside) const;

    Vec3 computeDifferenceWithCloestPoint(const Vec3& point, 
                                          const MatX& ppoint_pz,
                                          const int obs_id,
                                          MatX& pdiff_pz,
                                          double& isInside) const;

    Vec3 computeDifferenceWithCloestPoint(const Vec3& point, 
                                          const MatX& ppoint_pz,
                                          const Eigen::Array<MatX, 3, 1>& ppoint_pz_pz,
                                          const int obs_id,
                                          MatX& pdiff_pz,
                                          Eigen::Array<MatX, 3, 1>& pdiff_pz_pz,
                                          double& isInside) const;

    void computeDistance(const Vec3& point) override;

    void computeDistance(const Vec3& point, 
                         const MatX& ppoint_pz) override;

    void computeDistance(const Vec3& point, 
                         const MatX& ppoint_pz,
                         const Eigen::Array<MatX, 3, 1>& ppoint_pz_pz) override;

    // class members:
        // box information
    std::vector<Vec3> boxCenters;
    std::vector<Vec3> boxOrientation;
    std::vector<Vec3> boxSize;
    Eigen::Array<Mat3, 1, Eigen::Dynamic> boxR;

        // hyperplane representation
    Eigen::Array<Vec3, Eigen::Dynamic, Eigen::Dynamic> normals; //!< Normals of the box faces
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> intercepts; //!< Intercepts of the box faces

        // vertices representation
    Eigen::Array<Vec3, Eigen::Dynamic, Eigen::Dynamic> vertices; //!< Vertices of the box
};

}; // namespace RAPTOR

#endif // BOX_COLLISION_AVOIDANCE_H
