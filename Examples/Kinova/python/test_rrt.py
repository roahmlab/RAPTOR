import numpy as np
import matplotlib.pyplot as plt
import pybullet as p
import time
import scipy.io as sio

from kinova_dynamics import set_position

import sys
sys.path.append('/workspaces/RAPTOR/build/lib') # be careful about the path
import KinovaHLP_nanobind

### initializations
display_info = True
urdf_filename = "/workspaces/RAPTOR/Robots/kinova-gen3/gen3_2f85_fixed.urdf"

p.connect(p.GUI)
robot = p.loadURDF(urdf_filename, useFixedBase=True)
num_joints = p.getNumJoints(robot)

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
 
### high level planner (RRT) information
rrtplanner = KinovaHLP_nanobind.WaypointPlanningPybindWrapper(urdf_filename)

rrt_collision_buffer = 0.05 # m
include_gripper_or_not = True
rrtplanner_time = 5.0 # sec, time limit for RRT

rrtplanner.set_obstacles(obstacles, rrt_collision_buffer)
rrtplanner.set_start_goal(start, goal)
try:
    rrtpath = rrtplanner.plan(rrtplanner_time, include_gripper_or_not)
    print(len(rrtpath))
except Exception as e:
    print(e)
    rrtpath = []
    
# visualize RRT path in pybullet
for pos in rrtpath:
    set_position(robot, pos)
    time.sleep(0.1)
input("Press Enter to continue...")

p.disconnect()