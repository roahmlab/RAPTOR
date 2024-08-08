import pybullet as p
import pybullet_data as pd
import matplotlib.pyplot as plt
import time
import numpy as np

### settings
timeStep = 0.1

### read data
data = np.loadtxt("../data/trajectory-digit-simulation.txt")

### connect to simulator
p.connect(p.GUI)
p.setAdditionalSearchPath(pd.getDataPath())

# Load a simple plane
# plane_id = p.loadURDF("plane.urdf")
urdf_filename = "../../../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf"
robot = p.loadURDF(urdf_filename, useFixedBase=True)

# Start the simulation
p.setGravity(0, 0, -9.81)
p.setTimeStep(timeStep)
num_joints = p.getNumJoints(robot)

# input("Press Enter to continue...")

for tid in range(0, data.shape[1]):
    # base_xyz = data[0:3, tid]
    # base_rpy = data[3:6, tid]
    # base_quat = p.getQuaternionFromEuler(base_rpy)
    # pos = data[6:36, tid]
    
    # p.resetBasePositionAndOrientation(robot, base_xyz, base_quat)
    pos = data[:36, tid]
    
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

# Disconnect from PyBullet
p.disconnect()
