#ifndef BOX_SWEEP_VOLUME_COLLISION_AVOIDANCE_H
#define BOX_SWEEP_VOLUME_COLLISION_AVOIDANCE_H

#include "BoxCollisionAvoidance.h"

namespace IDTO {

namespace BoxSweepVolumeAndBox {
    constexpr int GENERATOR_NUM = 6;
    constexpr int HYPERPLANE_NUM = GENERATOR_NUM * (GENERATOR_NUM - 1);
};

/*
Check collision between the sweep volume of a box and other boxes.
*/
class BoxSweepVolumeCollisionAvoidance : public CollisionAvoidance {
public:
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    const int GENERATOR_NUM = BoxSweepVolumeAndBox::GENERATOR_NUM;
    const int HYPERPLANE_NUM = BoxSweepVolumeAndBox::HYPERPLANE_NUM;

    // Constructor
    BoxSweepVolumeCollisionAvoidance() = default;

    BoxSweepVolumeCollisionAvoidance(const Vec3& targetBoxCenter_input,
                                     const Vec3& targetBoxOrientation_input,
                                     const Vec3& targetBoxSize_input,
                                     const std::vector<Vec3>& boxCenters_input,
                                     const std::vector<Vec3>& boxOrientation_input,
                                     const std::vector<Vec3>& boxSize_input);

    // Destructor
    ~BoxSweepVolumeCollisionAvoidance() = default;

    // class methods:
    void initialize();

    void computeDistance(const Vec3& point) override;

    void computeDistance(const Vec3& point, const MatX& ppoint_pz) override;

    // class members:
        // target box information
    Vec3 targetBoxCenter;
    Vec3 targetBoxOrientation;
    Vec3 targetBoxSize;

        // other box information
    std::vector<Vec3> boxCenters;
    std::vector<Vec3> boxOrientation;
    std::vector<Vec3> boxSize;
    Eigen::Array<Mat3, 1, Eigen::Dynamic> boxR;

        // hyperplane representation
    Eigen::Array<Vec3, Eigen::Dynamic, Eigen::Dynamic> normals; //!< Normals of the box faces
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> intercepts; //!< Intercepts of the box faces

        // intermediate variables
    Eigen::VectorXd nominators;
    Eigen::VectorXd denominators;
    Eigen::VectorXd lambda_lb;
    Eigen::VectorXd lambda_ub;
};

}; // namespace IDTO

#endif // BOX_SWEEP_VOLUME_COLLISION_AVOIDANCE_H
