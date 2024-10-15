import pybullet as p
import pybullet_data
import numpy as np
import time
import sys
print(sys.executable)


def main():
    physicsClient = p.connect(p.GUI)
    p.setAdditionalSearchPath(pybullet_data.getDataPath())
    p.setGravity(0, 0, -9.81)

    planeId = p.loadURDF("plane.urdf")
    robotId = p.loadURDF("/home/RAPTOR/Robots/kinova-gen3/kinova.urdf", useFixedBase=1)
    numJoints = p.getNumJoints(robotId)

    with open("/home/RAPTOR/Examples/Kinova/SystemIdentification/position.csv", "r") as file:
        lines = file.readlines()

    trajectory = []

    for i, line in enumerate(lines):
        columns = line.split()  
        if len(columns) == 7:

            trajectory.append([float(x) for x in columns])
        else:
            print(f"Line {i+1} has {len(columns)} columns: {line}")

 
    print(f"Total trajectory steps: {len(trajectory)}")
        
    sphereRadius = 0.02  #
    sphereColor = [1, 0, 0, 1]  
    visualShapeId = p.createVisualShape(shapeType=p.GEOM_SPHERE, radius=sphereRadius, rgbaColor=sphereColor)

    for q_step in trajectory:
        for j in range(numJoints - 1):  # end is fixed
            p.resetJointState(robotId, j, q_step[j])
        p.stepSimulation()
        eeState = p.getLinkState(robotId, numJoints - 2)
        current_position = eeState[0]

        mass = 0  
        p.createMultiBody(baseMass=mass, baseVisualShapeIndex=visualShapeId, basePosition=current_position)

        print("End effector position:", current_position)
        time.sleep(1)

    print("Simulation finished, press Enter to exit...")
    input()
    p.disconnect()

if __name__ == "__main__":
    main()
