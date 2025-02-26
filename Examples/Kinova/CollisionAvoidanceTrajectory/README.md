# Discrete Version of ARMOUR

## Introduction
This folder contains an implementation of an discrete version of [ARMOUR](https://roahmlab.github.io/armour/).
Instead of guaranteeing safety on continuous time intervals, we only validate constraints on discrete time nodes.
This could be used to provide good initial guesses for [ARMOUR](https://roahmlab.github.io/armour/) because this runs much faster than [ARMOUR](https://roahmlab.github.io/armour/).

## Implementation Details
The following constraints are activated:
- Joint limits
- Velocity limits
- Torque limits
- Obstacle avoidance, where robot are represented as collection of spheres and the environment are represented as collection of boxes

Analytical Hessian of all constraints are provided.

The obstacles are defined in `KinovaExample.cpp`.
Check line 29-33.