# RAPTOR: RAPid and Robust Trajectory Optimization for Robots

## Introduction

Dynamic locomotion for humanoid robots presents significant analytical and computational challenges due to the extensive number of linkages and degrees of freedom. 
This complexity results in a vast search space for feasible gaits which translates into a time-consuming process when optimizing over trajectories. 
In addition, the process often involves numerous hyperparameters and requires a good initial guess or a warm-start strategy, further complicating the development process. 
Existing methods struggle to integrate the latest hardware designs, such as actuated ankles with closed-loop mechanisms, which offer increased stability, but introduce additional constraints into the dynamics that can be challenging to represent in a computationally tractable fashion. 
This work introduces a generalized gait optimization framework that directly generates smooth and physically feasible trajectories. 
The proposed method demonstrates faster and more robust convergence than existing techniques and explicitly incorporates closed-loop constraints. 
The method is implemented as an open-source C++ codebase that significantly reduces computation times, facilitating dynamic locomotion for full-size humanoids.
Note that **RAPTOR** also works for other fully actuated systems such as robotic manipulators.

To be more specific, we parameterize the trajectories for the robot states as, for example, polynomials.
The decision variable of the optimization is then the coefficients of the polynomial.
We sample a certain number of discrete points on the trajectory and evaluate the constraints, such as joints limits, torque limits, or collision avoidance.
We use Ipopt as our optimization solver.
For some of the constraints, we implement the analytical hessian so that Ipopt can converge faster.

This branch attempts to implement the sparse gradient for Ipopt. 
But it is unfinished right now and is not verified to be useful (significantly faster than dense gradient).

## Demos
<img src="https://github.com/user-attachments/assets/6f0a94cd-9c90-4d8f-ad6a-e7de86b017b6" alt="Digit walking forward" width="300" height="300">

<img src="https://github.com/user-attachments/assets/7c715902-3192-43ca-83a2-33239c758bf9" alt="Digit stepping stones" width="300" height="300">

<img src="https://github.com/user-attachments/assets/a68ab768-917e-4f6f-bd9b-64f0f337c025" alt="Talos walking forward" width="300" height="300">

## Requirements
- Ubuntu >= 20.04
- [Eigen 3.4](https://eigen.tuxfamily.org/index.php?title=3.4)
- [GSL](https://www.gnu.org/software/gsl/): solve close-loop kinematics
- [Pinocchio](https://stack-of-tasks.github.io/pinocchio/download.html): compute inverse dynamics and its gradient
- [Ipopt](https://coin-or.github.io/Ipopt/INSTALL.html): for nonlinear optimization

A more detailed instruction is provided [here](Installation/README.md).
We recommend users to install the requirements through docker.

## Overview
 - Trajectories/ : This folder contains implementation of multiple primitives of smooth trajectories.
 - KinematicsDynamics/ : This folder contains implementation to compute forward kinematics and inverse dynamics of a robot.
 - Constraints/ : This folder contains implementation of multiple constriants that could be useful for trajectory optimization,
                  such as torque limits or collision avoidance.
 - Optimization/ : This folder contains a base class that provides interfaces to ipopt.  
 - Examples/ : This folder contains two examples of trajectory optimization problem implementation for Digit-v3 and Kinova-gen3.

A more detailed instruction on how to code your own optimization problem is provided [here](Coding/README.md).  
            
## Getting Started
Run the following command to compile all the code.
```bash
mkdir build
cd build
cmake ..
make -j
```

You will be able to find README in each of the examples presented in `Example/` folder.
We provide the following examples
 - Kinova-gen3
    - `InverseKinematics/`: A simple inverse kinematics example given a desired end effector transformation matrix.
    - `ArmourDiscrete/`: The robot arm has nothing on its end effector. Reaching a target configuration while avoiding obstacles and satisfying torque limits.
    - `WaitrDiscrete/`: The robot arm is holding a tray with an object on it. Reaching a target configuratio while avoiding obstacles, satisfying torque limits, and making sure the object not fall off from the tray.
    - `SystemIdentification/`: **(Future work)** Multiple examples related to system identification of a robotic arm. 
 - Digit
    - Single step periodic gait optimization
    - Multiple step gait optimization
 - Talos
    - Single step gait optimization while starting from a fixed initial configuration
    - Multiple step gait optimization

## Bibtex
To cite **RAPTOR** in your academic research, please use the following bibtex entry:
```
@misc{zhang2024rapidrobusttrajectoryoptimization,
title={Rapid and Robust Trajectory Optimization for Humanoids}, 
author={Bohao Zhang and Ram Vasudevan},
year={2024},
eprint={2409.00303},
archivePrefix={arXiv},
primaryClass={cs.RO},
url={https://arxiv.org/abs/2409.00303}}
```

## Credits
Bohao Zhang (jimzhang@umich.edu)

This work is developed in [RoahmLab](http://www.roahmlab.com/), University of Michigan, Ann Arbor
