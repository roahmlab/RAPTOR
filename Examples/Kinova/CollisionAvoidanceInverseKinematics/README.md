# Kinova Inverse Kinematics Example

## Introduction

`KinovaIKExample.cpp` provides the simplest example of solving an inverse kinematics problem on Kinova.
Both analytical gradient and hessian are provided to Ipopt.

`KinovaIKMotionExample.cpp` provides an example of solving a series of inverse kinematics problems,
so that the end effector goes forward for 10 cm along its z axis.

## Implementation Details

The cost function is to minimize the distance from the initial guess.

The following constraints are activated:
- Joint limits
- Obstacle avoidance, where robot are represented as collection of spheres and the environment are represented as collection of boxes
- Kinematics constraints on the end effector position and orientation

The Hessian of this optimization problem is provided.