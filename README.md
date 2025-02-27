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

## Demos
<img src="https://github.com/user-attachments/assets/6f0a94cd-9c90-4d8f-ad6a-e7de86b017b6" alt="Digit walking forward" width="300" height="300">

<img src="https://github.com/user-attachments/assets/7c715902-3192-43ca-83a2-33239c758bf9" alt="Digit stepping stones" width="300" height="300">

<img src="https://github.com/user-attachments/assets/a68ab768-917e-4f6f-bd9b-64f0f337c025" alt="Talos walking forward" width="300" height="300">

<img src="https://github.com/user-attachments/assets/07b0763b-6042-49a2-94f5-a655f0553ad2" alt="Unitree G1 walking forward" width="300" height="300">

## Requirements
- Ubuntu >= 22.04
- [Eigen 3.4](https://eigen.tuxfamily.org/index.php?title=3.4)
- [GSL](https://www.gnu.org/software/gsl/): solve close-loop kinematics
- [Pinocchio >= 3.0](https://stack-of-tasks.github.io/pinocchio/download.html): compute inverse dynamics and its gradient
- [Ipopt](https://coin-or.github.io/Ipopt/INSTALL.html): for nonlinear optimization

A more detailed instruction is provided [here](Installation/README.md).
We recommend users to install the requirements through docker.

## Overview
 - Trajectories/ : This folder contains implementation of multiple primitives of smooth trajectories.
 - KinematicsDynamics/ : This folder contains implementation to compute forward kinematics and inverse dynamics of a robot.
 - Constraints/ : This folder contains implementation of multiple constriants that could be useful for trajectory optimization,
                  such as torque limits or collision avoidance.
 - Costs/ : This folder contains implementation of multiple costs that could be useful for trajectory optimization,
                  such as minimizing the total torque/power consumption or minimizing the path length.
 - Optimization/ : This folder contains a base class that provides interfaces to ipopt.  
 - Examples/ : This folder contains several examples of trajectory optimization problem implementations, including gait optimization examples for multiple walking robots like Digit-v3, Talos, Unitree-G1, and a lot of examples related to collision avoidance or system identification for Kinova-gen3.

A more detailed instruction on how to code your own optimization problem is provided [here](Coding/README.md).  
            
## Getting Started
Run the following command to compile all the code.
```bash
mkdir build
cd build
cmake ..
make -j4
```

You can find a README file with more details in each example within the `Example/` folder.
We provide the following examples:

- **Kinova-gen3**
    - `CollisionAvoidanceInverseKinematics/`: Solves an inverse kinematics problem for a given desired end-effector transformation matrix while considering joint limits and collision avoidance.
    - `CollisionAvoidanceTrajectory/`: The robot arm, without any payload, plans a trajectory to reach a target configuration while avoiding obstacles and satisfying torque limits.
    - `Armour`: A better re-implementation that integrates our previous work [ARMOUR](https://roahmlab.github.io/armour/) and [WAITR](https://roahmlab.github.io/waitr-dev/), which is a provably-safe receding-horizon trajectory optimization framework.
    - `SystemIdentification/`: Contains multiple examples related to system identification of a robotic arm.
    <!-- - `WaitrDiscrete/`: The robot arm holds a tray with an object on it. It optimizes a trajectory to reach a target configuration while avoiding obstacles, satisfying torque limits, and ensuring the object remains stable on the tray. -->

- **Digit**
    - Single-step periodic gait optimization.
    - Multi-step periodic gait optimization.
    - Detailed explanations on incorporating closed-loop linkages as kinematic constraints in optimization.

- **Talos**
    - Single-step periodic gait optimization.
    - Single-step gait optimization starting from a fixed initial configuration.
    - Multi-step periodic gait optimization.
    - Serves as a template for migrating the optimization framework to other humanoid robots.
    - Follow the procedure in the [README](Examples/Talos/README.md) of the folder to optimize gaits for your own humanoid robot.

- **Unitree-G1**
    - Single-step periodic gait optimization
    - Multi-step periodic gait optimization.
    - A simple example migrated from code in [Talos](Examples/Talos/).

## Bibtex
To cite **RAPTOR** in your academic research, please use the following bibtex entry:
```bibtex
@article{zhang2024rapidrobusttrajectoryoptimization,
  title={Rapid and Robust Trajectory Optimization for Humanoids},
  author={Bohao Zhang and Ram Vasudevan},
  journal={arXiv preprint arXiv:2409.00303},
  year={2024}
}
```

## Related Projects
**RAPTOR** also includes a better integrated implementation of our previous projects:

[Autonomous Robust Manipulation via Optimization with Uncertainty-aware Reachability](https://github.com/roahmlab/armour)
<!-- ```bibtex
@article{article,
author = {Michaux, Jonathan and Holmes, Patrick and Zhang, Bohao and Chen, Che and Wang, Baiyue and Sahgal, Shrey and Zhang, Tiancheng and Dey, Sidhartha and Kousik, Shreyas and Vasudevan, Ram},
year = {2023},
month = {01},
pages = {},
title = {Can't Touch This: Real-Time, Safe Motion Planning and Control for Manipulators Under Uncertainty}
doi={10.48550/arXiv.2301.13308}}
``` -->

[Wrench Analysis for Inertial Transport using Reachability](https://roahmlab.github.io/waitr-dev/)
<!-- ```bibtex
@ARTICLE{10403905,
  author={Brei, Zachary and Michaux, Jonathan and Zhang, Bohao and Holmes, Patrick and Vasudevan, Ram},
  journal={IEEE Robotics and Automation Letters}, 
  title={Serving Time: Real-Time, Safe Motion Planning and Control for Manipulation of Unsecured Objects}, 
  year={2024},
  volume={9},
  number={3},
  pages={2383-2390},
  keywords={Trajectory;Manipulators;Real-time systems;Uncertainty;Planning;Optimization;Safety;Collision avoidance;Manipulation planning;robot safety;collision avoidance},
  doi={10.1109/LRA.2024.3355731}}
``` -->

[Safe Planning for Articulated Robots Using Reachability-based Obstacle Avoidance With Spheres](https://roahmlab.github.io/sparrows/)
<!-- ```bibtex
@INPROCEEDINGS{Michaux-SPARROWS-RSS-24, 
  AUTHOR    = {Jonathan Michaux AND Adam Li AND Qingyi Chen AND Che Chen AND Ram Vasudevan},
  TITLE     = {{Safe Planning for Articulated Robots Using Reachability-based Obstacle Avoidance With Spheres}},
  BOOKTITLE = {Proceedings of Robotics: Science and Systems},
  YEAR      = {2024},
  ADDRESS   = {Delft, Netherlands},
  MONTH     = {July},
  DOI       = {10.15607/RSS.2024.XX.035}
}
``` -->

## Credits
[Bohao Zhang](https://cfather.github.io/) (jimzhang@umich.edu): Project leader.

Zichang Zhou (zhouzichang1234@gmail.com): System identification.

<!-- [Matthew Ejakov](https://sites.google.com/umich.edu/matthew-ejakov/home) (matthewejakov@gmail.com): Tapered capsule distance -->

Jiyang Wang (realwjy@umich.edu): Linearized reachability-based contact constraints. 

This work is developed in [RoahmLab](http://www.roahmlab.com/), University of Michigan, Ann Arbor
