import numpy as np
import matplotlib.pyplot as plt
import pybullet as p
import time
import scipy.io as sio

from kinova_dynamics import set_position

import sys
sys.path.append('/workspaces/RAPTOR/build/lib') # be careful about the path
import kinova_longer_nanobind

### initializations
display_info = True
urdf_filename = "/workspaces/RAPTOR/Robots/kinova-gen3/gen3_2f85_fixed.urdf"

p.connect(p.GUI)
robot = p.loadURDF(urdf_filename, useFixedBase=True)
num_joints = p.getNumJoints(robot)

# initialize planner
planner = kinova_longer_nanobind.KinovaLongerHorizonPybindWrapper(
    urdf_filename, 
    display_info
)

# obstacle information (xyz, rpy, size)
obstacles = np.array([[0.4, 0.0, 0.15, 0.0, 0.0, 0.0, 0.30, 0.30, 0.30], # another obstacle
                      [0.4, 0.0, 0.45, 0.0, 0.0, 0.0, 0.30, 0.30, 0.30], # another obstacle
                      [0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.01], # ground
                      [-0.39, -0.2, 0.56, 0.0, 0.0, 0.0, 0.08, 2.0, 1.12], # side wall
                     ])

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

# start and goal in configuration space
start = np.array([-2*np.pi/3, -1.0472, 0., -2.0944, 0., np.pi/2, 0.])
goal = np.array([2*np.pi/3, -1.0472, 0., -2.0944, 0., np.pi/2, 0.])

# initialize the robot to start configuration
set_position(robot, start)
 
# ### high level planner (RRT) information
# rrt_collision_buffer = 0.05 # m
# include_gripper_or_not = True
# rrtplanner_time = 5.0 # sec
# rrtplanner.set_obstacles(obstacles, rrt_collision_buffer)
# rrtplanner.set_start_goal(start, goal)
# try:
#     rrtpath = rrtplanner.plan(rrtplanner_time, include_gripper_or_not)
#     print(len(rrtpath))
# except Exception as e:
#     print(e)
#     rrtpath = []
    
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
# set up local planner (dump all the information into it)
planner.set_ipopt_parameters(
    1e-5,          # ipopt_tol
    1e-6,          # ipopt_constr_viol_tol
    10.0,          # ipopt_obj_scaling_factor
    20.0,          # ipopt_max_wall_time
    5,             # ipopt_print_level
    "adaptive",    # mu_strategy
    "ma57",        # ipopt_linear_solver   
    False          # ipopt_gradient_check
)

collision_buffer = 0.02 # m
planner.set_obstacles(
    obstacles,
    collision_buffer
)

# trajectory information
duration = 4.0 # sec, duration of the trajectory
degree = 2 # number of control points in the middle
planner.set_trajectory_parameters(
    start, 
    goal, 
    duration, 
    degree
)

# buffer information (to make sure constraints are satisfied between time nodes)
joint_limits_buffer = np.array([0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02])
velocity_limits_buffer = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
torque_limits_buffer = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
planner.set_buffer(
    joint_limits_buffer, \
    velocity_limits_buffer, \
    torque_limits_buffer
)

z, if_success = planner.optimize()

if if_success:
    traj, _, _, _ = planner.analyze_solution()

    # replay the trajectories in pybullet
    for tid in range(traj.shape[0]):
        set_position(robot, traj[tid, :7])
        time.sleep(0.05)
        
    input("Press Enter to continue...")
else:
    raise ValueError("Optimization failed")

p.disconnect()