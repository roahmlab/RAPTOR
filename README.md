# IDTO
Trajectory Optimization based on Inverse Dynamics

## Requirements
- Ubuntu >= 20.04
- [Eigen 3.4](https://eigen.tuxfamily.org/index.php?title=3.4)
- [Pinocchio](https://stack-of-tasks.github.io/pinocchio/download.html)
- [Ipopt](https://coin-or.github.io/Ipopt/INSTALL.html)
- PkgConfig (`pkg-config`)
- `urdfdom` (`liburdfdom-dev`)
- [DRD](https://github.com/Cfather/DRD)

A more detailed instruction is provided in the `Installation/README.md` file.

## Overview
 - Trajectories/ : This folder contains implementation of multiple primitives of smooth trajectories.
 - InverseDynamics/ : This folder contains implementation to compute forward kinematics and inverse dynamics of a robot.
 - Constraints/ : This folder contains implementation of multiple constriants that could be useful for trajectory optimization,
                  such as torque limits or collision avoidance.
 - Optimization/ : This folder contains a base class that provides interfaces to ipopt.                 

## Getting Started
Run the following commands for a planar 3DOF robot example:
```bash
mkdir build
cd build
cmake ..
make
./IDTO_example
```

## Credits
Bohao Zhang (jimzhang@umich.edu)

RoahmLab, University of Michigan, Ann Arbor