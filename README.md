# IDTO
Trajectory Optimization based on Inverse Dynamics

We aim to solve a trajectory optimization problem formulated as below

![ProblemFormulation](Assets/pic-ProblemFormulation.svg)

where we assume that the trajectories are already parameterized and only treat the trajectory parameters as the decision variables.

For example, you can define your trajectory to be a polynomial (on a fixed interval [0,T]):

![TrajectoryFormulation](Assets/pic-TrajectoryFormulation.svg)

And based on that, you can add joint limit constraints, torque limit constraints, end effector constraints, or any other customized constraints to your optimization problem.

These constraints will then be evaluated on a sequence of discrete time instances t_i over [0,T], so that the trajectory respects all of the constraints over all time.

Compared with the direct collocation method or the single/multiple shooting method, this method guarantees a **smooth** trajectory instead of a discrete representation.
And this method could usually run faster since less decision variables are involved.

Note that this formulation only works for fully-actuated systems.

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
 - Examples/ : This folder contains two examples of trajectory optimization problem for Digit-v3 and Kinova-gen3.                

## Getting Started
Run the following commands for a Digit-v3 example:
```bash
mkdir build
cd build
cmake ..
make
./IDTO_example
```

## Credits
Bohao Zhang (jimzhang@umich.edu)

[RoahmLab](http://www.roahmlab.com/), University of Michigan, Ann Arbor