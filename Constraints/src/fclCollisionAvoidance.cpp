#include "fclCollisionAvoidance.h"

namespace IDTO {

bool CustomizedDistanceFunction(fcl::CollisionObjectd* o1, 
                                fcl::CollisionObjectd* o2, 
                                void* cdata_, 
                                double& dist) {
    assert(cdata_ != nullptr);
    auto* cdata = static_cast<CustomizedDistanceData*>(cdata_);
    const fcl::DistanceRequestd& request = cdata->request;
    fcl::DistanceResultd& result = cdata->result;

    if (cdata->done) { 
        dist = result.min_distance; 
        return true; 
    }

    const fcl::NODE_TYPE o1type = o1->getNodeType(); 
    const fcl::NODE_TYPE o2type = o2->getNodeType();
    dist = 0;
    std::string name1, name2;

    if (o1type == fcl::GEOM_BOX && 
        o2type == fcl::GEOM_SPHERE) {
        const auto sphereGeomtry = static_cast<const fcl::Sphered*>(o2->collisionGeometry().get());
        const double sphere_radius = sphereGeomtry->radius;

        dist = computeBoxPointDistance(o1, o2->getTranslation()) - sphere_radius;

        const auto userData1 = (customizedUserDataForBox*)(o1->getUserData());
        const auto userData2 = (customizedUserDataForSphere*)(o2->getUserData());
        name1 = userData1->name;
        name2 = userData2->name;
    }
    else if (o1type == fcl::GEOM_SPHERE && 
             o2type == fcl::GEOM_BOX) {
        const auto sphereGeomtry = static_cast<const fcl::Sphered*>(o1->collisionGeometry().get());
        const double sphere_radius = sphereGeomtry->radius;

        dist = computeBoxPointDistance(o2, o1->getTranslation()) - sphere_radius;

        const auto userData1 = (customizedUserDataForSphere*)(o1->getUserData());
        const auto userData2 = (customizedUserDataForBox*)(o2->getUserData());
        name1 = userData1->name;
        name2 = userData2->name;
    }
    else if (o1type == fcl::GEOM_SPHERE && 
             o2type == fcl::GEOM_SPHERE) {
        const auto userData1 = (customizedUserDataForSphere*)(o1->getUserData());
        const auto userData2 = (customizedUserDataForSphere*)(o2->getUserData());

        // we don't check collision between adajacent links
        if (std::abs(userData1->linkId - userData2->linkId) <= 1) {
            dist = std::numeric_limits<double>::max();
        }
        else {
            const auto sphereGeomtry1 = static_cast<const fcl::Sphered*>(o1->collisionGeometry().get());
            const auto sphereGeomtry2 = static_cast<const fcl::Sphered*>(o2->collisionGeometry().get());
            const double sphere_radius1 = sphereGeomtry1->radius;
            const double sphere_radius2 = sphereGeomtry2->radius;

            dist = (o1->getTranslation() - o2->getTranslation()).norm() - (sphere_radius1 + sphere_radius2);
        }

        name1 = userData1->name;
        name2 = userData2->name;
    }
    else {
        std::cerr << "The collision objects should be box-sphere or sphere-sphere!\n";
        std::cerr << "Other differentiable collision checking is not supported yet.\n";
        throw std::invalid_argument("The collision objects should be box-sphere or sphere-sphere!");
    }

    if (dist < result.min_distance) {
        result.min_distance = dist;
        cdata->min_distance = dist;
        cdata->name1 = name1;
        cdata->name2 = name2;
    }

    if (dist <= 0) {
        return true; // in collision or in touch
    }

    return cdata->done;
}

bool CustomizedDistanceFunctionDerivative(fcl::CollisionObjectd* o1, 
                                          fcl::CollisionObjectd* o2, 
                                          void* cdata_, 
                                          double& dist) {
    assert(cdata_ != nullptr);
    auto* cdata = static_cast<CustomizedDistanceData*>(cdata_);
    const fcl::DistanceRequestd& request = cdata->request;
    fcl::DistanceResultd& result = cdata->result;

    if (cdata->done) { 
        dist = result.min_distance; 
        return true; 
    }

    const fcl::NODE_TYPE o1type = o1->getNodeType(); 
    const fcl::NODE_TYPE o2type = o2->getNodeType();
    dist = 0;
    Eigen::VectorXd pdistance_pz;
    std::string name1, name2;

    if (o1type == fcl::GEOM_BOX && 
        o2type == fcl::GEOM_SPHERE) {
        const auto userData1 = (customizedUserDataForBox*)(o1->getUserData());
        const auto userData2 = (customizedUserDataForSphere*)(o2->getUserData());
        const Eigen::MatrixXd& ppoint_pz = userData2->ppoint_pz;
        const auto sphereGeomtry = static_cast<const fcl::Sphered*>(o2->collisionGeometry().get());
        const double sphere_radius = sphereGeomtry->radius;

        std::pair<double, Eigen::VectorXd> distRes = 
            computeBoxPointDistance(o1, o2->getTranslation(), ppoint_pz);

        dist = distRes.first - sphere_radius;
        pdistance_pz = distRes.second;

        name1 = userData1->name;
        name2 = userData2->name;
    }
    else if (o1type == fcl::GEOM_SPHERE && 
             o2type == fcl::GEOM_BOX) {
        const auto userData1 = (customizedUserDataForSphere*)(o1->getUserData());
        const auto userData2 = (customizedUserDataForBox*)(o2->getUserData());
        const Eigen::MatrixXd& ppoint_pz = userData1->ppoint_pz;
        const auto sphereGeomtry = static_cast<const fcl::Sphered*>(o1->collisionGeometry().get());
        const double sphere_radius = sphereGeomtry->radius;

        std::pair<double, Eigen::VectorXd> distRes = 
            computeBoxPointDistance(o2, o1->getTranslation(), ppoint_pz);
        
        dist = distRes.first - sphere_radius;
        pdistance_pz = distRes.second;
        
        name1 = userData1->name;
        name2 = userData2->name;
    }
    else if (o1type == fcl::GEOM_SPHERE && 
             o2type == fcl::GEOM_SPHERE) {
        const auto userData1 = (customizedUserDataForSphere*)(o1->getUserData());
        const auto userData2 = (customizedUserDataForSphere*)(o2->getUserData());

        // we don't check collision between adajacent links
        if (std::abs(userData1->linkId - userData2->linkId) <= 1) {
            dist = std::numeric_limits<double>::max();
        }
        else {
            const Eigen::MatrixXd& ppoint_pz1 = userData1->ppoint_pz;
            const Eigen::MatrixXd& ppoint_pz2 = userData2->ppoint_pz;
            const auto sphereGeomtry1 = static_cast<const fcl::Sphered*>(o1->collisionGeometry().get());
            const auto sphereGeomtry2 = static_cast<const fcl::Sphered*>(o2->collisionGeometry().get());
            const double sphere_radius1 = sphereGeomtry1->radius;
            const double sphere_radius2 = sphereGeomtry2->radius;

            dist = (o1->getTranslation() - o2->getTranslation()).norm();

            if (dist > 1e-5) {
                pdistance_pz = (o1->getTranslation() - o2->getTranslation()).transpose() * 
                                    (ppoint_pz1 - ppoint_pz2) / dist;
            }
            else {
                pdistance_pz.setZero();
            }

            dist -= (sphere_radius1 + sphere_radius2);
        }

        name1 = userData1->name;
        name2 = userData2->name;
    }
    else {
        std::cerr << "The collision objects should be box-sphere or sphere-sphere!\n";
        std::cerr << "Other differentiable collision checking is not supported yet.\n";
        throw std::invalid_argument("The collision objects should be box-sphere or sphere-sphere!");
    }

    if (dist < result.min_distance) {
        result.min_distance = dist;
        cdata->min_distance = dist;
        cdata->name1 = name1;
        cdata->name2 = name2;
        cdata->pmin_distance_pz = pdistance_pz;
    }

    if (dist <= 0) {
        return true; // in collision or in touch
    }

    return cdata->done;
}

void fclCollisionAvoidance::addObstacleBox(const std::string& name, 
                                           const Vec3& boxCenter,
                                           const Vec3& boxOrientation,
                                           const Vec3& boxSize) {
    std::shared_ptr<fcl::Boxd> boxGeomtry(new fcl::Boxd(boxSize(0), 
                                                        boxSize(1), 
                                                        boxSize(2)));

    Mat3 R = (Eigen::AngleAxisd(boxOrientation[0], Vec3::UnitX())
                * Eigen::AngleAxisd(boxOrientation[1], Vec3::UnitY())
                * Eigen::AngleAxisd(boxOrientation[2], Vec3::UnitZ())).matrix();
    
    fclObjectCollection[name] = new fcl::CollisionObjectd(boxGeomtry, 
                                                          R, 
                                                          boxCenter);

    fclObjectCollection[name]->setUserData((void*)(new customizedUserDataForBox(name)));
    fclObjectCollection[name]->computeAABB();

    fclBroadPhaseManagerPtr_->registerObject(fclObjectCollection[name]);
    fclBroadPhaseManagerPtr_->setup();
}

void fclCollisionAvoidance::addRobotSphere(const std::string& name, 
                                           const int linkId,
                                           const Vec3& sphereCenter,
                                           const double sphereRadius,
                                           const MatX& ppoint_pz) {
    std::shared_ptr<fcl::Sphered> sphereGeomtry(new fcl::Sphered(sphereRadius));
    fclObjectCollection[name] = new fcl::CollisionObjectd(sphereGeomtry,
                                                          Mat3::Identity(),
                                                          sphereCenter);

    fclObjectCollection[name]->setUserData((void*)(new customizedUserDataForSphere(name, linkId, ppoint_pz)));
    fclObjectCollection[name]->computeAABB();

    fclBroadPhaseManagerPtr_->registerObject(fclObjectCollection[name]);
    fclBroadPhaseManagerPtr_->setup();
}

void fclCollisionAvoidance::clear() {
    fclBroadPhaseManagerPtr_->clear();

    for (const auto& object : fclObjectCollection) {
        if (object.second->getNodeType() == fcl::NODE_TYPE::GEOM_BOX) {
            delete (customizedUserDataForBox*)(object.second->getUserData());
        }
        else if (object.second->getNodeType() == fcl::NODE_TYPE::GEOM_SPHERE) {
            delete (customizedUserDataForSphere*)(object.second->getUserData());
        }
        else {
            throw std::runtime_error("The collision object should be a box or a sphere!");
        }
        delete object.second;
    }

    fclObjectCollection.clear();
}

void fclCollisionAvoidance::computeDistance(bool compute_derivatives) {
    distanceData.reset();
    if (compute_derivatives) {
        fclBroadPhaseManagerPtr_->distance(&distanceData, 
                                           &CustomizedDistanceFunctionDerivative);
    }
    else {
        fclBroadPhaseManagerPtr_->distance(&distanceData, 
                                           &CustomizedDistanceFunction);
    }
}

void fclCollisionAvoidance::computeDistance(const std::shared_ptr<fclCollisionAvoidance>& otherCollisionManager,
                                            bool compute_derivatives) {
    distanceData.reset();
    if (compute_derivatives) {
        fclBroadPhaseManagerPtr_->distance(otherCollisionManager->fclBroadPhaseManagerPtr_.get(), 
                                           &distanceData, 
                                           &CustomizedDistanceFunctionDerivative);
    }
    else {
        fclBroadPhaseManagerPtr_->distance(otherCollisionManager->fclBroadPhaseManagerPtr_.get(), 
                                           &distanceData, 
                                           &CustomizedDistanceFunction);
    }
}

}; // namespace IDTO