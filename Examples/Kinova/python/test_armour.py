import numpy as np
import sys
import matplotlib.pyplot as plt
import time
import pybullet as p

from kinova_dynamics import goal_distance, set_position

sys.path.append("/workspaces/RAPTOR/build/lib")
import armour_nanobind as armour
import KinovaHLP_nanobind as HLP

### initializations
display_info = False

urdf_filename = "/workspaces/RAPTOR/Robots/kinova-gen3/gen3_2f85_fixed.urdf" # this urdf only has 7 joints for visualization
config_filename = "/workspaces/RAPTOR/Examples/Kinova/Armour/KinovaWithGripperInfo.yaml" # edit this yaml file if necessary

p.connect(p.GUI)
robot = p.loadURDF(urdf_filename, useFixedBase=True)
num_joints = p.getNumJoints(robot)

rrtplanner = HLP.WaypointPlanningPybindWrapper(urdf_filename)
planner = armour.ArmourPybindWrapper(urdf_filename, config_filename, display_info)

# obstacle information (xyz, rpy, size)
obstacles = np.array([[0.3, 0, 0.25, 0.0, 0.0, 0.1, 0.10, 0.10, 0.50], # obstacle in the middle
                      [0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.01], # ground
                      [0.53, 0.49, 0.56, 0.0, 0.0, 0.0, 2.0, 0.08, 1.12], # back wall
                      [-0.39, -0.82, 0.56, 0.0, 0.0, 0.0, 0.08, 0.08, 1.12], # camera bar near the control
                      [-0.39, 0.42, 0.56, 0.0, 0.0, 0.0, 0.08, 0.08, 1.12], # second camera bar
                      [0.7, 0.0, 1.12, 0.0, 0.0, 0.0, 2.0, 2.0, 0.05] # ceiling
                     ])

# start and goal in joint configuration space
start = np.array([1.28159142, 0.8, -2.6481 , -2.3484, 2.55565072, -0.83405794, 2.05711487])
goal = np.array([-0.8, 0.4, -2.8913, -2.0785, 2.86698255, 0.41219441, -0.2324])

# initialize the robot to start configuration
set_position(robot, start)

# create and visualize the obstacles in pybullet
for i, obs in enumerate(obstacles):
    position = obs[0:3]
    orientation = p.getQuaternionFromEuler(obs[3:6])
    size = obs[6:9]
    
    boxCollisionShapeId = p.createCollisionShape(shapeType=p.GEOM_BOX,
                                                 halfExtents=size*0.5)
    boxVisualShapeId = p.createVisualShape(shapeType=p.GEOM_BOX,
                                           halfExtents=size*0.5,
                                           rgbaColor=[1,0,0,0.2])
    box_id = p.createMultiBody(baseMass=0,
                               baseVisualShapeIndex=boxVisualShapeId,
                               baseCollisionShapeIndex=boxCollisionShapeId,
                               basePosition=position,
                               baseOrientation=orientation)
 
 
# high level planner (RRT) information
rrt_collision_buffer = 0.02 # m
include_gripper_or_not = True
rrtplanner_time = 5.0 # sec, timeout for RRT planner
rrtplanner.set_obstacles(obstacles, rrt_collision_buffer)
rrtplanner.set_start_goal(start, goal)
try:
    rrtpath = rrtplanner.plan(rrtplanner_time, include_gripper_or_not)
    print(len(rrtpath))
except Exception as e:
    print(e)
    rrtpath = []
    
# # visualize RRT path in pybullet
# for pos in rrtpath:
#     id = 0
#     for i in range(num_joints):
#         joint_info = p.getJointInfo(robot, i)
#         joint_type = joint_info[2]
#         if joint_type == p.JOINT_FIXED:
#             p.resetJointState(robot, i, targetValue=0)
#         else:
#             p.resetJointState(robot, i, targetValue=pos[id])
#             id += 1
#     p.stepSimulation()
#     time.sleep(0.1)
# input("Press Enter to continue...")

### start local planning
duration = 2.0 # sec, duration of the (Bezier) trajectory

# set up local planner (dump all the information into it)
planner.set_obstacles(obstacles)
planner.set_ipopt_parameters(
    1e-6,           # tol
    1e-6,           # constr_viol_tol
    10.0,           # obj_scaling_factor
    5.0,            # max_wall_time
    5,              # print_level
    "adaptive",     # mu_strategy
    "ma57",         # linear_solver
    False           # gradient_check
)    

# initial conditions of the robot   
q0 = start
qd0 = np.zeros(7)
qdd0 = np.zeros(7)
 
# choose to minimize the distance between the goal and which part of the trajectory
# for example, 0.0 means minimizing the distance between the goal and the start of the trajectory
# 1.0 means minimizing the distance between the goal and the end of the trajectory
t_plan = 1.0 # choose between [0.0, 1.0]

for iter in range(20):
    # set up the local planner
    q_des = goal # desired joint configuration of the current planning iteration
    k_center = np.pi/12 * (q_des - q0) / np.linalg.norm(q_des - q0) # searching k at this point, pi/6 forward from start to goal
    k_range = np.pi/24 * np.ones(7)

    planner.set_trajectory_parameters(
        q0, \
        qd0, \
        qdd0, \
        k_center, \
        k_range, \
        duration, \
        q_des, \
        t_plan
    )

    # optimize the trajectory
    k, ifFeasible = planner.optimize()
    print(k)

    if ifFeasible:
        info, \
        sphere_x, sphere_y, sphere_z, sphere_radii, \
        torque_centers, torque_radii = \
            planner.analyze_solution()
        
        sphere_info = np.stack((sphere_x, sphere_y, sphere_z), axis=0)
        previous_info = info
            
        # visualize the trajectory
        for tid in range(info.shape[1]):
            set_position(robot, info[:7, tid])
            time.sleep(0.01)
    else:
        print("Failed to find a feasible trajectory")
        break
    
    # update robot information for the next planning iteration
    q0 = info[:7, -1]
    qd0 = np.zeros(7)
    qdd0 = np.zeros(7)
    
    if goal_distance(q0 - goal) < 0.1:
        print("Goal reached")
        break
        
p.disconnect()