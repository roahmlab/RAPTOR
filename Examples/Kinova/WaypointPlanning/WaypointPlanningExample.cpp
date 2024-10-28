#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/base/ScopedState.h>
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/ProblemDefinition.h>
#include <ompl/geometric/SimpleSetup.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/config.h>

#include "KinovaConstants.h"
#include "Plain.h"
#include "KinovaCustomizedConstraints.h"

#include <iostream>

// OMPL namespaces for easier use
namespace ob = ompl::base;
namespace og = ompl::geometric;

using namespace RAPTOR;
using namespace Kinova;

int main() {
    std::srand(std::time(nullptr));

    // Initialize the robot model
    const std::string urdf_filename = "../Robots/kinova-gen3/kinova.urdf";
    
    pinocchio::Model model_double;
    pinocchio::urdf::buildModel(urdf_filename, model_double);
    pinocchio::ModelTpl<double> model = model_double.cast<double>();

    // Initialize obstacles
    const int num_obstacles = 12;
    std::vector<Eigen::Vector3d> boxCenters;
    std::vector<Eigen::Vector3d> boxOrientation;
    std::vector<Eigen::Vector3d> boxSize;

    boxCenters.resize(num_obstacles);
    boxOrientation.resize(num_obstacles);
    boxSize.resize(num_obstacles);

    boxCenters[0] << -0.20, 1.25, 0.85;
    boxCenters[1] << -0.20, -1.25, 0.85; 
    boxCenters[2] << -0.35, 0.25, 0.85; 
    boxCenters[3] << -0.35, -0.25, 0.85; 
    boxCenters[4] << -0.35, 0.00, 1.10; 
    boxCenters[5] << -0.35, 0.00, 0.60; 
    boxCenters[6] << -0.50, 0.00, 0.85; 
    boxCenters[7] << 0.60, -0.30, 0.85; 
    boxCenters[8] << 0.60, -0.80, 0.85; 
    boxCenters[9] << 0.85, -0.55, 0.85; 
    boxCenters[10] << 0.35, -0.55, 0.85; 
    boxCenters[11] << 0.60, -0.55, 1.10; 

    for (int i = 0; i < num_obstacles; i++) {
        boxOrientation[i].setZero();
    }

    boxSize[0] << 0.01, 1.00, 0.25;
    boxSize[1] << 0.01, 1.00, 0.25; 
    boxSize[2] << 0.15, 0.01, 0.25; 
    boxSize[3] << 0.15, 0.01, 0.25; 
    boxSize[4] << 0.15, 0.25, 0.01; 
    boxSize[5] << 0.15, 0.25, 0.01; 
    boxSize[6] << 0.01, 0.25, 0.25; 
    boxSize[7] << 0.25, 0.01, 0.25; 
    boxSize[8] << 0.25, 0.01, 0.25; 
    boxSize[9] << 0.01, 0.25, 0.25; 
    boxSize[10] << 0.01, 0.25, 0.25; 
    boxSize[11] << 0.25, 0.25, 0.01; 

    // Initialize the collision checker
    std::shared_ptr<Trajectories> trajPtr_ = std::make_shared<Plain>(NUM_JOINTS);
    std::shared_ptr<Constraints> collisionCheckerPtr_ = 
        std::make_shared<KinovaCustomizedConstraints>(
            trajPtr_,
            model,
            boxCenters,
            boxOrientation,
            boxSize
        );

    // Construct the robot configuration space
    auto space = std::make_shared<ob::RealVectorStateSpace>(NUM_JOINTS);

    // Set the bounds of the space
    ob::RealVectorBounds bounds(NUM_JOINTS);
    for (int i = 0; i < NUM_JOINTS; i++) {
        const double lower_bound = 
            (JOINT_LIMITS_LOWER[i] == -1e19) ? 
                -M_PI : 
                JOINT_LIMITS_LOWER[i];

        const double upper_bound = 
            (JOINT_LIMITS_UPPER[i] == 1e19) ? 
                M_PI : 
                JOINT_LIMITS_UPPER[i];

        bounds.setLow(i, lower_bound);
        bounds.setHigh(i, upper_bound);
    }
    space->setBounds(bounds);

    // Create a simple setup object
    og::SimpleSetup ss(space);

    // Set state validity checker
    const double buffer = 0.0;
    ss.setStateValidityChecker(
        [&collisionCheckerPtr_, &buffer](const ob::State *state) { 
            Eigen::VectorXd joint_angles(NUM_JOINTS);
            for (int i = 0; i < NUM_JOINTS; i++) {
                joint_angles(i) = state->as<ob::RealVectorStateSpace::StateType>()->values[i];
            }

            collisionCheckerPtr_->compute(joint_angles, false);
            
            for (int i = 0; i < collisionCheckerPtr_->m; i++) {
                if (collisionCheckerPtr_->g(i) <= buffer) {
                    return false;
                }
            }

            return true;
        });

    // Define the start and goal states
    ob::ScopedState<> start(space);
    start[0] = 0;
    start[1] = 0.5236;
    start[2] = 0;
    start[3] = -1.1972;
    start[4] = 0;
    start[5] = -1.0472;
    start[6] = 0;

    ob::ScopedState<> goal(space);
    goal[0] = 0.5236;
    goal[1] = 1.309;
    goal[2] = -1.5708;
    goal[3] = -0.3927;
    goal[4] = 1.5708;
    goal[5] = -1.5708;
    goal[6] = 0;

    ss.setStartAndGoalStates(start, goal);

    // Use RRTConnect as the planner
    auto planner = std::make_shared<og::RRTConnect>(ss.getSpaceInformation());
    planner->setRange(0.05);

    ss.setPlanner(planner);

    // Attempt to solve the problem within 1 second
    ob::PlannerStatus solved = ss.solve(1.0);

    if (solved) {
        std::cout << "Found solution:" << std::endl;
        // ss.simplifySolution();
        ss.getSolutionPath().print(std::cout);
    }
    else {
        std::cout << "No solution found" << std::endl;
    }

    return 0;
}
