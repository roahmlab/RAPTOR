# Discrete Version of WAITR

This folder contains an implementation of an discrete version of [WAITR](https://roahmlab.github.io/waitr-dev/).
Instead of guaranteeing safety on continuous time intervals, we only validate constraints on discrete time nodes.
This could be used to provide good initial guesses for [WAITR](https://roahmlab.github.io/waitr-dev/) because this runs much faster than [WAITR](https://roahmlab.github.io/waitr-dev/).

The following constraints are activated:
- Joint limits
- Velocity limits
- Torque limits
- Obstacle avoidance, where robot are represented as collection of spheres and the environment are represented as collection of boxes
- Contact constraints (circle contact surface) between the tray and the object

The obstacles are defined in `KinovaWaitrExample.cpp`.
Check line 36-40.