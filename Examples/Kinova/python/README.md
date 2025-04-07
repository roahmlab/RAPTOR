# Python Scripts for Kinova

This folder contains python scripts related to multiple trajectory optimization examples and system identification examples for Kinova.

## Before You Run These Scripts

### Path to run the scripts
Make sure you are inside the current folder (`RAPTOR/Examples/Kinova/python`) to run these scripts because of the path defined in these scripts.

### Path to the URDF and the pybind program
We assume that you are running everything inside the docker.
That's why in these scripts, you will find that all the paths are defined globally inside the docker, starting with `/workspaces/RAPTOR/`.
If you run the scripts without the docker, you will need to manually change the paths to all the URDFs and the pybind programs.

### Other Notes
In these scripts, we use PyBullet **only for visualization**.
The dynamics of the robot is simulated separately using `scipy.integrate.solve_ivp` for solving ODE and `pinocchio.aba` for forward dynamics.
This approach is generally much more accurate than physics simulators using Euler integration method, like PyBullet.
Since collision shapes of the obstacles are enabled, you should make sure that the robot is not in collision with any obstacles in PyBullet.
Otherwise, the visualization of the robot might be incorrect due to the collision.

## Introduction to Each Script

### kinova_dynamics.py
A collection of helper functions for the Kinova robot that provides:
- Robot model loading and initialization
- Forward and inverse dynamics computations
- Integration utilities for simulation and visualization in PyBullet

### visualize_trajectory.py
A basic example that animates the robot trajectories inside PyBullet.
In line 141-167 of `KinovaExample.cpp`, you will find that we usually export the trajectory as a text file in the following format:
```C++
std::ofstream trajectory("trajectory-kinova.txt");
trajectory << std::setprecision(20);
for (int i = 0; i < NUM_JOINTS; i++) {
    for (int j = 0; j < N; j++) {
        trajectory << mynlp->trajPtr_->q(j)(i) << ' ';
    }
    trajectory << std::endl;
}
for (int i = 0; i < NUM_JOINTS; i++) {
    for (int j = 0; j < N; j++) {
        trajectory << mynlp->trajPtr_->q_d(j)(i) << ' ';
    }
    trajectory << std::endl;
}
for (int i = 0; i < NUM_JOINTS; i++) {
    for (int j = 0; j < N; j++) {
        trajectory << mynlp->trajPtr_->q_dd(j)(i) << ' ';
    }
    trajectory << std::endl;
}
for (int i = 0; i < NUM_JOINTS; i++) {
    for (int j = 0; j < N; j++) {
        trajectory << mynlp->idPtr_->tau(j)(i) << ' ';
    }
    trajectory << std::endl;
}
trajectory.close();
```
As a result, you can read the text file in `visualize_trajectory.py` and animate the motion by replaying the joint positions in sequence.

line 17-39 `visualize_trajectory.py` also shows how to defined box obstacles and visualize in PyBullet, which is applied in all example scripts in this folder.
Each box obstacle is represented by a 9-element array with the following structure:
```python
obstacles = np.array([
    [x, y, z, roll, pitch, yaw, width, height, depth], # 1st obstacle
    [x, y, z, roll, pitch, yaw, width, height, depth], # 2nd obstacle
    ...
    [x, y, z, roll, pitch, yaw, width, height, depth], # nth obstacle
])
```
where:
- `[x, y, z]`: Position of the box center in meters
- `[roll, pitch, yaw]`: Orientation of the box in radians (converted to quaternion internally)
- `[width, height, depth]`: The full dimensions of the box in meters

The obstacles are visualized in PyBullet with:
- Semi-transparent red color (rgba=[1,0,0,0.2])
- Zero mass and fixed (static obstacles not considered in dynamics)
- Both collision and visual shapes for proper physics simulation and visualization

### test_sysid_momentum.py
System identification script that:
- Collect robot trajectory data from a dynamics simulation of robot tracking a sine wave using a naive PD controller.
- Momentum-based system identification for the Kinova robot based on the collect data. Note that this method **does not require robot acceleration**.
- Outputs estimates of end-effector inertial parameters and **the uncertainties of the estimates**

### test_sysid_inverse_dynamics.py
Similar to `test_sysid_momentum.py` but uses inverse dynamics approach for system identification.
Note that this method **requires robot acceleration**.

### test_rrt.py
Interface to our implementations of an RRT planner based on [ompl](https://ompl.kavrakilab.org/).

### test_kinova_longer_horizon_planner.py
An example script of generating safe and smooth trajectories, that satisfy the following constraints over **a predefined set of discrete time instances**:

- Joint limits
- Velocity limits
- Torque limits
- Collision avoidance with box obstacles

This is mostly common used to generate global collision avoidance trajectories in complex scenarios, while the constraints might be violated between the discrete time instances.

More information can be found at the [README](../CollisionAvoidanceTrajectory/README.md) in `CollisionAvoidanceTrajectory` folder.

### test_armour.py
An example script to interface our previous work [ARMOUR](https://roahmlab.github.io/armour/).
[ARMOUR](https://roahmlab.github.io/armour/) can generate safe and smooth trajectories in a receding horizon fashion, that satisfy the following constraints over **continuous time intervals**:

- Joint limits
- Velocity limits
- Torque limits
- Collision avoidance with box obstacles

More information can be found at the [README](../Armour/README.md) in `Armour` folder.