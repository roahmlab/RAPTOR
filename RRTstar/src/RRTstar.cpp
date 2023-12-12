#include "RRTstar.h"

namespace IDTO {

RRTStar::RRTStar(const VecX& start_input, 
                 const VecX& goal_input, 
                 const VecX& lowerLimits_input,
                 const VecX& upperLimits_input, 
                 const RRTParameters& params_input,
                 std::shared_ptr<PlanarBoxObstacleAvoidance>& collisionPtr_input) :
    start(start_input), 
    goal(goal_input), 
    lowerLimits(lowerLimits_input),
    upperLimits(upperLimits_input), 
    params(params_input),
    collisionPtr_(collisionPtr_input) {
    if (start.size() != goal.size()) {
        throw std::invalid_argument("Start and goal dimensions do not match!");
    }

    if (start.size() != lowerLimits.size()) {
        throw std::invalid_argument("Start and lower limits dimensions do not match!");
    }

    if (start.size() != upperLimits.size()) {
        throw std::invalid_argument("Start and upper limits dimensions do not match!");
    }

    if (params.stepSize <= 0) {
        throw std::invalid_argument("Step size must be positive!");
    }

    if (params.maxIterations <= 0) {
        throw std::invalid_argument("Maximum number of iterations must be positive!");
    }

    // Initialize the tree with the start node
    nodes_.push_back(std::make_shared<Node>(start));

    trajPtr_ = std::make_shared<Plain>(start.size());
    collisionPtr_->trajPtr_ = trajPtr_;

    collisionPtr_->compute_bounds();
}

RRTStar::~RRTStar() {
    deleteNodes();
}

void RRTStar::deleteNodes() {
    // for (Node*& node : nodes_) {
    //     delete node;
    //     node = nullptr; // Set pointer to null after deletion
    // }
    // nodes_.clear();
}

std::vector<Eigen::VectorXd> RRTStar::plan() {
    for (int i = 0; i < params.maxIterations; ++i) {
        if (i % 1000 == 0) {
            int report = nearestNode(goal);
            std::cout << "Iter: " << i << " Current distance to goal: " << distance(nodes_[report]->point, goal) << std::endl;
        }

        // Generate a random point
        VecX randomPoint = generateRandomNode();

        // Find the nearest node in the tree
        int nearest_id = nearestNode(randomPoint);
        std::shared_ptr<Node>& nearest = nodes_[nearest_id];

        // Steer towards the random point
        VecX newPoint = steer(nearest->point, randomPoint);

        // Check if the path is collision-free
        if (isCollisionFree(nearest->point, newPoint)) {
            // Create a new node
            std::shared_ptr<Node> newNode = std::make_shared<Node>(newPoint);
            newNode->parent = nearest;

            // std::cout << "New node: " << newNode->point.transpose() << std::endl;

            // Find nearby nodes
            // std::vector<Node*> nearNodes;
            // for (Node* node : nodes_) {
            //     if (distance(node->point, newPoint) < 2 * stepSize) {
            //         nearNodes.push_back(node);
            //     }
            // }

            // Connect the new node to the tree
            nearest->children.push_back(newNode);

            // Rewire the tree
            // rewire(newNode, nearNodes);

            // Add the new node to the list
            nodes_.push_back(newNode);

            // Check if the goal is reached
            if (distance(newPoint, goal) <= 2 * params.stepSize) {
                std::cout << "Goal reached!" << std::endl;

                // Construct the path from the goal to the start
                std::vector<VecX> path;
                std::shared_ptr<Node> current = newNode;
                while (current != nullptr) {
                    path.push_back(current->point);
                    current = current->parent;
                }

                // Reverse the path to get it from start to goal
                std::reverse(path.begin(), path.end());

                return path;
            }
        }
    }

    // If the goal is not reached, return an empty path
    return std::vector<VecX>();
}

Eigen::VectorXd RRTStar::generateRandomNode() {
    // Still choose goal with with a small probability
    if (static_cast<double>(rand()) / RAND_MAX < params.forwardGoalProbability) {
        return goal;
    }

    // Generate a random point within the workspace
    VecX randomPoint(start.size());
    for (int i = 0; i < start.size(); ++i) {
        randomPoint(i) = lowerLimits(i) + (static_cast<double>(rand()) / RAND_MAX) * (upperLimits(i) - lowerLimits(i));
    }

    return randomPoint;
}

int RRTStar::nearestNode(const VecX& randomPoint) {
    // Find the nearest node in the tree
    int nearest_id = -1;
    double minDistance = std::numeric_limits<double>::max();

    for (int i = 0; i < nodes_.size(); i++) {
        double d = distance(nodes_[i]->point, randomPoint);
        if (d < minDistance) {
            nearest_id = i;
            minDistance = d;
        }
    }

    return nearest_id;
}

Eigen::VectorXd RRTStar::steer(const VecX& from, const VecX& to) {
    // Steer towards the target point
    double d = distance(from, to);
    if (d < params.stepSize) {
        return to;
    } 
    else {
        VecX direction = (to - from).normalized();
        return from + params.stepSize * direction;
    }
}

bool RRTStar::isCollisionFree(const VecX& from, const VecX& to) {
    // Check if the line segment from 'from' to 'to' is collision-free
    
    // Only check from and to since we assume we use a very small step size
    collisionPtr_->compute(from, false);
    for (int i = 0; i < collisionPtr_->m; i++) {
        if (collisionPtr_->g(i) < collisionPtr_->g_lb(i) || 
            collisionPtr_->g(i) > collisionPtr_->g_ub(i)) {
            return false;
        }
    }

    collisionPtr_->compute(to, false);
    for (int i = 0; i < collisionPtr_->m; i++) {
        if (collisionPtr_->g(i) < collisionPtr_->g_lb(i) || 
            collisionPtr_->g(i) > collisionPtr_->g_ub(i)) {
            return false;
        }
    }

    return true;
}

void RRTStar::rewire(Node* newNode, const std::vector<std::shared_ptr<Node>>& nearNodes) {
    // Rewire the tree by checking if the new node provides a shorter path to nearby nodes
    // for (int i = 0; i < nearNodes.size(); i++) {
    //     const std::shared_ptr<Node>& nearNode = nearNodes[i];
    //     double cost = distance(newNode->point, nearNode->point) + nearNode->point[0];

    //     if (cost < distance(nearNode->point, newNode->parent->point) + newNode->parent->point[0]) {
    //         // Rewire the tree
    //         newNode->parent->children.erase(
    //             std::remove(newNode->parent->children.begin(), newNode->parent->children.end(), nearNode),
    //             newNode->parent->children.end());
    //         nearNode->parent = newNode;
    //         newNode->parent->children.push_back(nearNode);
    //     }
    // }
}

double RRTStar::distance(const VecX& p1, const VecX& p2) {
    // Euclidean distance between two points
    return (p1 - p2).norm();
}

}; // namespace IDTO
