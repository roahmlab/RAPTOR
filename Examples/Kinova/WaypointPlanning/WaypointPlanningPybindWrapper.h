#ifndef KINOVA_HLP_PYBIND_WRAPPER_H
#define KINOVA_HLP_PYBIND_WRAPPER_H

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

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

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

namespace RAPTOR {
namespace Kinova {

namespace nb = nanobind;
namespace ob = ompl::base;
namespace og = ompl::geometric;

class WaypointPlanningPybindWrapper {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    using nb_1d_double = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
    using nb_2d_double = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

    // Constructor
    WaypointPlanningPybindWrapper() = default;

    WaypointPlanningPybindWrapper(const std::string urdf_filename) {
        pinocchio::Model model_double;
        pinocchio::urdf::buildModel(urdf_filename, model_double);
        model = model_double.cast<double>();

        // Construct the robot configuration space
        space = std::make_shared<ob::RealVectorStateSpace>(NUM_JOINTS);

        // Set the bounds of the space
        ob::RealVectorBounds bounds(NUM_JOINTS);
        for (int i = 0; i < NUM_JOINTS; i++) {
            const double lower_bound = 
                (JOINT_LIMITS_LOWER[i] == -1e19) ? 
                    -4 * M_PI : 
                    JOINT_LIMITS_LOWER[i] * 0.9;

            const double upper_bound = 
                (JOINT_LIMITS_UPPER[i] == 1e19) ? 
                    4 * M_PI : 
                    JOINT_LIMITS_UPPER[i] * 0.9;

            bounds.setLow(i, lower_bound);
            bounds.setHigh(i, upper_bound);
        }
        space->setBounds(bounds);

        // Create a simple setup object
        ss = std::make_shared<og::SimpleSetup>(space);

        auto planner = std::make_shared<og::RRTConnect>(ss->getSpaceInformation());
        planner->setRange(0.05);

        ss->setPlanner(planner);
    };

    // Destructor
    ~WaypointPlanningPybindWrapper() = default;

    // Class methods
    void set_obstacles(const nb_2d_double obstacles_inp,
                       const double buffer_inp) {
        if (obstacles_inp.shape(1) != 9) {
            throw std::invalid_argument("Obstacles must have 9 columns, xyz, rpy, size");
        }

        if (buffer_inp < 0) {
            throw std::invalid_argument("Collision checking buffer must be non-negative");
        }

        buffer = buffer_inp;

        num_obstacles = obstacles_inp.shape(0);

        boxCenters.resize(num_obstacles);
        boxOrientation.resize(num_obstacles);
        boxSize.resize(num_obstacles);

        for (int i = 0; i < num_obstacles; i++) {
            boxCenters[i] << obstacles_inp(i, 0), obstacles_inp(i, 1), obstacles_inp(i, 2);
            boxOrientation[i] << obstacles_inp(i, 3), obstacles_inp(i, 4), obstacles_inp(i, 5);
            boxSize[i] << obstacles_inp(i, 6), obstacles_inp(i, 7), obstacles_inp(i, 8);
        }

        set_obstacles_check = true;
    };

    void set_start_goal(const nb_1d_double start_inp, 
                        const nb_1d_double goal_inp) {
        if (start_inp.shape(0) != NUM_JOINTS || goal_inp.shape(0) != NUM_JOINTS) {
            throw std::invalid_argument("Start and goal must be of size NUM_JOINTS");
        }

        for (int i = 0; i < NUM_JOINTS; i++) {
            if (start_inp(i) < JOINT_LIMITS_LOWER[i] || start_inp(i) > JOINT_LIMITS_UPPER[i]) {
                std::cerr << i << ' ' << start_inp(i) << " [" << JOINT_LIMITS_LOWER[i] << ' ' << JOINT_LIMITS_UPPER[i] << "]\n"; 
                throw std::invalid_argument("Start state is out of joint limits");
            }

            if (goal_inp(i) < JOINT_LIMITS_LOWER[i] || goal_inp(i) > JOINT_LIMITS_UPPER[i]) {
                std::cerr << i << ' ' << goal_inp(i) << " [" << JOINT_LIMITS_LOWER[i] << ' ' << JOINT_LIMITS_UPPER[i] << "]\n";
                throw std::invalid_argument("Goal state is out of joint limits");
            }
        }

        // Define the start and goal states
        ob::ScopedState<> start(space);
        ob::ScopedState<> goal(space);

        for (int i = 0; i < NUM_JOINTS; i++) {
            start[i] = start_inp(i);
            goal[i] = goal_inp(i);
        }

        ss->setStartAndGoalStates(start, goal);

        set_start_goal_check = true;
    };

    nb::ndarray<nb::numpy, const double> plan(const double timeout,
                                              const bool include_gripper_or_not = false) {
        if (!set_obstacles_check ||
            !set_start_goal_check) {
            throw std::runtime_error("Obstacles and start/goal must be set before planning");
        }

        std::shared_ptr<Trajectories> trajPtr_ = 
            std::make_shared<Plain>(NUM_JOINTS);
        std::shared_ptr<KinovaCustomizedConstraints> collisionCheckerPtr_ =
            std::make_shared<KinovaCustomizedConstraints>(
                trajPtr_,
                model,
                boxCenters,
                boxOrientation,
                boxSize,
                include_gripper_or_not
            );

        const double buffer_local = buffer;

        ss->setStateValidityChecker(
        [&collisionCheckerPtr_, &buffer_local](const ob::State *state) { 
            Eigen::VectorXd joint_angles(NUM_JOINTS);
            for (int i = 0; i < NUM_JOINTS; i++) {
                joint_angles(i) = state->as<ob::RealVectorStateSpace::StateType>()->values[i];
            }

            collisionCheckerPtr_->compute(joint_angles, false);
            
            for (int i = 0; i < collisionCheckerPtr_->m; i++) {
                if (collisionCheckerPtr_->g(i) <= buffer_local) {
                    return false;
                }
            }

            return true;
        });

        ob::PlannerStatus solved = ss->solve(timeout);

        if (!solved) {
            throw std::runtime_error("Failed to find a solution");
        }

        // ss->simplifySolution();
        auto sol = ss->getSolutionPath();
        sol.interpolate();

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> path(sol.getStateCount(), NUM_JOINTS);
        for (int i = 0; i < sol.getStateCount(); i++) {
            auto state = sol.getState(i);
            for (int j = 0; j < NUM_JOINTS; j++) {
                path(i, j) = state->as<ob::RealVectorStateSpace::StateType>()->values[j];
            }
        }

        set_start_goal_check = false;

        const size_t shape_ptr[] = {sol.getStateCount(), NUM_JOINTS};
        return nb::ndarray<nb::numpy, const double>(
            path.data(),
            2,
            shape_ptr,
            nb::handle()
        );
    };

    // Class members
    // robot model
    Model model;

    // obstacle information
    int num_obstacles = 0;
    std::vector<Vec3> boxCenters;
    std::vector<Vec3> boxOrientation;
    std::vector<Vec3> boxSize;
    double buffer = 0.0;

    // planner information
    std::shared_ptr<ob::RealVectorStateSpace> space;
    std::shared_ptr<og::SimpleSetup> ss;

    // Flags to check if the parameters are set
    bool set_obstacles_check = false;
    bool set_start_goal_check = false;
};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // KINOVA_HLP_PYBIND_WRAPPER_H
