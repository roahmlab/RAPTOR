import numpy as np
import time
import sys
import matplotlib.pyplot as plt


with open("/home/zichang/RAPTOR/Examples/Kinova/SystemIdentification/ExcitingTrajectories/data/T10_d5_slower/exciting-velocity-11.csv", "r") as file:
    lines = file.readlines()

velocity_data = [
    [float(x) for x in line.split()]
    for line in lines if len(line.split()) == 7
]

with open("/home/zichang/RAPTOR/Examples/Kinova/SystemIdentification/ExcitingTrajectories/data/T10_d5_slower/exciting-acceleration-11.csv", "r") as file:
    lines = file.readlines()

acceleration_data = [
    [float(x) for x in line.split()]
    for line in lines if len(line.split()) == 7
]

num_joints = len(acceleration_data[0])
time_steps = range(len(acceleration_data))

fig, axes = plt.subplots(nrows=num_joints, ncols=2, figsize=(16, num_joints * 4))

# Plot data for all joints
for joint in range(num_joints):
    # Plot Velocity
    ax_vel = axes[joint, 0]
    ax_vel.plot(time_steps, [step[joint] for step in velocity_data], color='orange')
    ax_vel.set_ylabel("Velocity (rad/s)")
    ax_vel.set_title(f"Joint {joint + 1} Velocity")
    ax_vel.grid(True)

    # Plot Acceleration
    ax_acc = axes[joint, 1]
    ax_acc.plot(time_steps, [step[joint] for step in acceleration_data], color='red')
    ax_acc.set_xlabel("Time Steps")
    ax_acc.set_ylabel("Acceleration (rad/s²)")
    ax_acc.set_title(f"Joint {joint + 1} Acceleration")
    ax_acc.grid(True)

# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()