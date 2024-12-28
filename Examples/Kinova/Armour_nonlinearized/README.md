# ARMOUR
This folder contains an implementation of [ARMOUR](https://roahmlab.github.io/armour/).

The following constraints are activated:
- Joint limits
- Velocity limits
- Torque limits
- Obstacle avoidance, where robot are represented as collection of spheres and the environment are represented as collection of boxes

Analytical hessian of all constraints are provided.

<!-- The obstacles are defined in `ArmourExample.cpp`.
Check line 29-33. -->