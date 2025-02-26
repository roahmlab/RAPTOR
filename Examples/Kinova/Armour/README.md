# ARMOUR
This folder contains an implementation that integrates both [ARMOUR](https://roahmlab.github.io/armour/) and [WAITR](https://roahmlab.github.io/waitr-dev/).

## Introduction

A key challenge to the widespread deployment of robotic manipulators is the need to ensure safety in arbitrary environments while generating new motion plans in real-time. 
In particular, one must ensure that the manipulator does not collide with obstacles, collide with itself, or exceed its own joint torque limits. 
This challenge is compounded by the need to account for uncertainty in the mass and inertia of manipulated objects, and potentially the robot itself. 
The present work addresses this challenge by proposing Autonomous Robust Manipulation via Optimization with Uncertainty-aware Reachability **ARMOUR**, a provably-safe, receding-horizon trajectory planner and tracking controller framework for serial link manipulators. 

In particular, this paper makes three contributions. 
1. a robust, passivity-based controller enables a manipulator to track desired trajectories with bounded error despite uncertain dynamics. 
2. a novel variation on the Recursive Newton-Euler Algorithm (RNEA) allows **ARMOUR** to compute the set of all possible inputs required to track any trajectory within a continuum of desired trajectories. 
3. this paper provides a method to compute the swept volume of the manipulator given a reachable set of states; this enables one to guarantee safety by checking that the swept volume does not intersect with obstacles. 

The proposed method is compared to state-of-the-art methods and demonstrated on a variety of challenging manipulation examples in simulation, such as maneuvering a heavy dumbbell with uncertain mass around obstacles.

## Implementation Details

The following constraints are activated:
- Joint limits
- Velocity limits
- Torque limits
- Obstacle avoidance, where robot are represented as collection of spheres and the environment are represented as collection of boxes

The Hessian of this optimization problem is not provided.