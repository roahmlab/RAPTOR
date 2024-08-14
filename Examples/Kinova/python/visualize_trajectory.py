import pybullet as p
import pybullet_data as pd
import matplotlib.pyplot as plt
import time
import numpy as np

### settings
timeStep = 0.1

### read data
data = np.loadtxt("../../../trajectory.txt")

### connect to simulator
p.connect(p.GUI)
p.setAdditionalSearchPath(pd.getDataPath())

# Load a simple plane
# plane_id = p.loadURDF("plane.urdf")
urdf_filename = "../../../Robots/kinova-gen3/kinova_grasp.urdf"
robot = p.loadURDF(urdf_filename, useFixedBase=True)

# Start the simulation
p.setGravity(0, 0, -9.81)
p.setTimeStep(timeStep)
num_joints = p.getNumJoints(robot)

# input("Press Enter to continue...")

for tid in range(data.shape[0]):
    for i in range(data.shape[1]):
        p.resetJointState(robot, i, targetValue=data[tid, i])
    
    p.stepSimulation()
    time.sleep(1e-2)
    
input("Press Enter to continue...")

# Disconnect from PyBullet
p.disconnect()
