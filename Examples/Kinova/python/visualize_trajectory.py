import numpy as np
import pickle
import matplotlib.pyplot as plt
import pybullet as p
import time
import math
import scipy.io as sio

from kinova_dynamics import set_position

urdf_filename = "/workspaces/RAPTOR/Robots/kinova-gen3/gen3_2f85_fixed.urdf"

p.connect(p.GUI)
robot = p.loadURDF(urdf_filename, useFixedBase=True)
num_joints = p.getNumJoints(robot)

# obstacle information (xyz, rpy, size)
obstacles = np.array([[0.40, 0.40, 0.15, 0.0, 0.0, 0.0, 0.30, 0.30, 0.30], # another obstacle
                      [0.40, 0.40, 0.45, 0.0, 0.0, 0.0, 0.30, 0.30, 0.30], # another obstacle
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
 
# load the trajectory from the text file
trajectory = np.loadtxt("/workspaces/RAPTOR/Examples/Kinova/CollisionAvoidanceTrajectory/data/T4_d2_obs/trajectory.txt")

for tid in range(trajectory.shape[0]):
    set_robot_position(trajectory[tid, :num_joints])
    time.sleep(0.05)
    
input("Press Enter to continue...")
p.disconnect()