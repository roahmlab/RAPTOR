#ifndef RRTSTAR_H
#define RRTSTAR_H

#include "Plain.h" 
#include "PlanarBoxObstacleAvoidance.h"

#include <vector>
#include <limits>
#include <memory>
#include <cstdlib>

namespace IDTO {

typedef struct RRTParameters_ {
    double stepSize = 0.1;
    int maxIterations = 10000;
    double forwardGoalProbability = 0.3;
    double collisionCheckThreshold = 0.02; // meter
} RRTParameters;

typedef struct Node_ {
    Eigen::VectorXd point;
    // std::vector<Node_*> children;
    // Node_* parent;
    std::vector<std::shared_ptr<Node_>> children;
    std::shared_ptr<Node_> parent;

    Node_(const Eigen::VectorXd& point) :
        point(point), parent(nullptr) {}

    ~Node_() {
        // for (Node_* child : children) {
        //     delete child; // Recursively delete children
        //     child = nullptr; // Set pointer to null after deletion
        // }
    }
} Node;

class RRTStar {
public:
    using VecX = Eigen::VectorXd;

    // Constructor
    RRTStar() = default;

    // Constructor
    RRTStar(const VecX& start_input, 
            const VecX& goal_input, 
            const VecX& lowerLimits_input,
            const VecX& upperLimits_input, 
            const RRTParameters& params_input,
            std::shared_ptr<PlanarBoxObstacleAvoidance>& collisionPtr_input);

    // Destructor
    ~RRTStar();

    // class methods:
    void deleteNodes();

    std::vector<VecX> plan();

    VecX generateRandomNode();

    int nearestNode(const VecX& randomPoint);

    VecX steer(const VecX& from, const VecX& to);

    bool isCollisionFree(const VecX& from, const VecX& to);

    void rewire(Node* newNode, const std::vector<std::shared_ptr<Node>>& nearNodes);

    double distance(const VecX& p1, const VecX& p2);

    // class variables:
    VecX start;
    VecX goal;
    
    VecX lowerLimits;
    VecX upperLimits;

    RRTParameters params;

    std::vector<std::shared_ptr<Node>> nodes_;

    std::shared_ptr<Plain> trajPtr_;

    // Construct a collision checker instance outside of the RRTStar class
    std::shared_ptr<PlanarBoxObstacleAvoidance> collisionPtr_;
};

}; // namespace IDTO

#endif // RRTSTAR_H