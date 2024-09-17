#ifndef FCL_COLLISION_AVOIDANCE_H
#define FCL_COLLISION_AVOIDANCE_H

#include <fcl/broadphase/default_broadphase_callbacks.h>
#include <fcl/narrowphase/collision_object.h>
#include <fcl/geometry/shape/sphere.h>
#include <fcl/narrowphase/distance.h>
#include <fcl/broadphase/broadphase_dynamic_AABB_tree.h>

#include "BoxCollisionAvoidance.h"

namespace RAPTOR {

struct customizedUserDataForSphere {
    customizedUserDataForSphere() {}
    customizedUserDataForSphere(const std::string name_input,
                                const int link_id_input) : 
        name(name_input),
        linkId(link_id_input){
        ppoint_pz = Eigen::MatrixXf::Zero(0, 0);
    }
    customizedUserDataForSphere(const std::string name_input,
                                const int link_id_input,
                                const Eigen::MatrixXf& ppoint_pz_input) : 
        name(name_input),
        linkId(link_id_input),
        ppoint_pz(ppoint_pz_input) {}

    std::string name = ""; //!< Name of the object
    int linkId = 0; //!< The link that this sphere belongs to (we don't check collision between adajacent links)
    Eigen::MatrixXf ppoint_pz; //!< Derivative of the center of the sphere (this is only for sphere!)
};

struct customizedUserDataForBox {
    customizedUserDataForBox() {}
    customizedUserDataForBox(const std::string name_input) : 
        name(name_input) {
    }

    std::string name = ""; //!< Name of the object
};

struct CustomizedDistanceData {
    CustomizedDistanceData() { 
        result.min_distance = std::numeric_limits<float>::max();
        min_distance = std::numeric_limits<float>::max();
        done = false; 
    }

    void reset() {
        result.min_distance = std::numeric_limits<float>::max();
        min_distance = std::numeric_limits<float>::max();
        done = false;
    }

    // the following are default values
    fcl::DistanceRequestd request; //!< Request parameters for distance computation.
    fcl::DistanceResultd result;   //!< Result of distance computation.
    bool done = false;                            //!< Flag indicating if distance computation is completed.

    // the following are customized values
    std::string name1 = "";               //!< Name of the first object in the collision pair.
    std::string name2 = "";               //!< Name of the second object in the collision pair.
    float min_distance = 0.0;            //!< Minimum distance between two collision managers.
    Eigen::VectorXf pmin_distance_pz;     //!< Derivative of the minimum distances between two collision managers.
};

bool CustomizedDistanceFunction(fcl::CollisionObjectd* o1, 
                                fcl::CollisionObjectd* o2, 
                                void* cdata_, 
                                float& dist);

bool CustomizedDistanceFunctionDerivative(fcl::CollisionObjectd* o1, 
                                          fcl::CollisionObjectd* o2, 
                                          void* cdata_, 
                                          float& dist);

class fclCollisionAvoidance {
public:
    using Vec3 = Eigen::Vector3f;
    using Mat3 = Eigen::Matrix3f;
    using VecX = Eigen::VectorXf;
    using MatX = Eigen::MatrixXf;

    fclCollisionAvoidance() {
        fclBroadPhaseManagerPtr_ = std::make_shared<fcl::DynamicAABBTreeCollisionManagerd>();
    }

    ~fclCollisionAvoidance() {
        clear();
    }

    bool empty() const {
        return fclObjectCollection.empty();
    }

    void addObstacleBox(const std::string& name, 
                        const Vec3& boxCenter,
                        const Vec3& boxOrientation,
                        const Vec3& boxSize);

    void addRobotSphere(const std::string& name, 
                        const int linkId,
                        const Vec3& sphereCenter,
                        const float sphereRadius,
                        const MatX& ppoint_pz = MatX::Zero(0, 0));

    void clear();

    // self collision check
    void computeDistance(bool compute_derivatives = false);

    // collision check with other collision manager
    void computeDistance(const std::shared_ptr<fclCollisionAvoidance>& otherCollisionManager,
                         bool compute_derivatives = false);

    std::unordered_map<std::string, fcl::CollisionObjectd*> fclObjectCollection; //!< Collection of collision objects.

    std::shared_ptr<fcl::BroadPhaseCollisionManagerd> fclBroadPhaseManagerPtr_ = nullptr; //!< Broad phase collision manager.

    CustomizedDistanceData distanceData; //!< Data for customized distance computation.
};

}; // namespace RAPTOR

#endif // FCL_COLLISION_AVOIDANCE_H