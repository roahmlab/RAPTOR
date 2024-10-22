import pybullet as p
import pybullet_data
import numpy as np
import time
import sys
import matplotlib.pyplot as plt
import math


print(sys.executable)

def create_obstacles(centers, sizes):
    """Create visual objects in PyBullet based on centers and sizes."""
    obstacle_visuals = []
    for center, size in zip(centers, sizes):
        # Define half extents since PyBullet uses half sizes for box shapes
        half_extents = [dim / 2.0 for dim in size]
        visual_shape_id = p.createVisualShape(
            shapeType=p.GEOM_BOX, 
            halfExtents=half_extents, 
            rgbaColor=[0, 1, 0, 0.5]  # Green with transparency
        )
        # Create the obstacle at the specified position
        obstacle_id = p.createMultiBody(
            baseMass=0,  # No physics for visual-only objects
            baseVisualShapeIndex=visual_shape_id,
            basePosition=center
        )
        obstacle_visuals.append(obstacle_id)
    return obstacle_visuals

def main():




    physicsClient = p.connect(p.GUI)
    p.setAdditionalSearchPath(pybullet_data.getDataPath())
    p.setGravity(0, 0, -9.81)

    # Load the plane and the robot URDF
    planeId = p.loadURDF("plane.urdf")
    robotId = p.loadURDF("/home/zichang/RAPTOR/Robots/kinova-gen3/kinova.urdf", useFixedBase=1)

    numJoints = 7  # Kinova Gen3 typically has 7 joints

    # Load the trajectory data from the CSV
    # with open("/home/zichang/RAPTOR/Examples/Kinova/SystemIdentification/ExcitingTrajectories/data/T10_d5_slower/exciting-position-18.csv", "r") as file:
    with open("/home/zichang/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/q_downsampled_5.csv", "r") as file:

        lines = file.readlines()

    trajectory = [
        [float(x) for x in line.split()]
        for line in lines if len(line.split()) == 7
    ]
    print(f"Total trajectory steps: {len(trajectory)}")

    # Define centers and sizes based on your C++ code
    centers = [
        [0.0, 0.0, 0.15],  # Floor center
        [0.53, 0.49, 0.56],  # Back wall
        [-0.39, -0.84, 0.56],  # Bar near control
        [-0.39, -0.17, 0.56],  # Bar between 10 and 20 (wall)
        [0.0, 0.0, 1.12]  # Ceiling
    ]
    sizes = [
        [5.0, 5.0, 0.01],  # Floor
        [5.0, 0.08, 1.12],  # Back wall
        [0.08, 0.08, 1.12],  # Bar near control
        [0.10, 1.28, 1.28],  # Wall between 10 and 20
        [5.0, 5.0, 0.05]  # Ceiling
    ]

    # Create obstacles in the scene
    # create_obstacles(centers, sizes)


    # Visualize trajectory execution
    sphereRadius = 0.01
    sphereColor = [1, 0, 0, 1]  
    visualShapeId = p.createVisualShape(p.GEOM_SPHERE, radius=sphereRadius, rgbaColor=sphereColor)

    # q_step = trajectory[0]

    q_step = [139.139, 75.747, 192.08, 318.838, 2.543, 24.144, 91.86]
    q_step = [math.radians(angle) for angle in q_step] 
    for j in range(numJoints):
        p.resetJointState(robotId, j, q_step[j])
    p.stepSimulation()
    eeState = p.getLinkState(robotId, numJoints - 1)  # Last link state (end effector)
    current_position = eeState[0]
    print("End effector position:", current_position)


    # for q_step in trajectory:
    #     for j in range(numJoints):
    #         p.resetJointState(robotId, j, q_step[j])
    #     p.stepSimulation()
    #     eeState = p.getLinkState(robotId, numJoints - 1)  # Last link state (end effector)
    #     current_position = eeState[0]

    #     p.createMultiBody(baseMass=0, baseVisualShapeIndex=visualShapeId, basePosition=current_position)

    #     print("End effector position:", current_position)
  

    print("Simulation finished, press Enter to exit...")
    input()
    p.disconnect()

if __name__ == "__main__":
    main()
