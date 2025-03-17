import pybullet as p
import pybullet_data as pd
import matplotlib.pyplot as plt
import time
import numpy as np

# load the trajectory data here
data = np.loadtxt("../data/full-trajectories-upstairs.txt")
data = data.T

# connect to simulator
p.connect(p.GUI)
p.setAdditionalSearchPath(pd.getDataPath())

# create and visualize the robot in pybullet
urdf_filename = "../../../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf"
robot = p.loadURDF(urdf_filename, useFixedBase=True)
nq = 36 # number of joints

# start the simulation
p.setGravity(0, 0, -9.81)
num_joints = p.getNumJoints(robot)

# visualize the trajectory along the time
for tid in range(0, data.shape[1]):
    pos = data[:nq, tid]
    
    # directly set the joint angles
    id = 0
    for i in range(num_joints):
        joint_info = p.getJointInfo(robot, i)
        joint_type = joint_info[2]
        if joint_type == p.JOINT_FIXED:
            p.resetJointState(robot, i, targetValue=0)
        else:
            p.resetJointState(robot, i, targetValue=pos[id])
            id += 1
    
    p.stepSimulation()
    time.sleep(1e-2)
    
input("Press Enter to continue...")

p.disconnect()
