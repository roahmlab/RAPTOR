# Robot Model Parameters Identification

## Introduction
This folder holds multiple examples for identification of robot model parameters:

### Motor friction parameters identification

Identifying motor friction parameters is crucial for accurate dynamic modeling and control of robotic manipulators. Unlike the rigid links of the robot, whose inertial parameters remain constant over time, motors can degrade or develop defects. Consequently, motor friction may need to be identified frequently to maintain precise torque estimation, optimize control performance, and ensure smooth motion execution.

`TestFrictionParametersIdentification.cpp` provides an example of identifying robot motor friction parameters.
This assumes that the robot inertial parameters have been identified and updated in the URDF.
The following motor friction model is used:
```math
F(\dot{q}(t), \ddot{q}(t), \theta) = F_{c} \circ \text{sign}(\dot{q}(t)) + F_{v} \circ \dot{q}(t) + I_{a} \circ \ddot{q}(t) + \beta
```

 - $F_c$: static friction coefficients.
 - $F_v$: viscous friction coefficients (motor damping coefficients).
 - $I_a$: transmission inertias (motor armature coefficients).
 - $\beta$: offset/bias terms.

### End effector inertial parameters identification

Robotic arms are widely used in cooperative human-robot environments, including manufacturing, package delivery, and in-home care. These scenarios often involve manipulating payloads with uncertain properties while operating under physical constraints.

`TestEndEffectorParametersIdentification.cpp` provides an example of identifying the inertial parameters of the end effector, assuming that the inertial parameters of the remaining robot structure (all links except the end effector) have already been identified and updated in the URDF.

This script applies the Log-Cholesky parameterization technique, first introduced in [1], to ensure the physical consistency of the identified inertial parameters. As a result, the optimization problem is reformulated as a nonlinear regression, which can be efficiently solved using Ipopt.

### End effector inertial parameters identification based on momentum

`TestEndEffectorParametersIdentificationMomentum.cpp` provides an alternative approach to identifying the end-effector inertial parameters, leveraging the momentum dynamics of the system. For more details on momentum dynamics, refer to Section II in [2].

This formulation eliminates the need for acceleration estimation by integrating the rigid body dynamics, which is beneficial since joint acceleration estimation is often intractable in most robotic systems. As a result, this approach is generally more robust to sensor noise compared to the traditional method used in `TestEndEffectorParametersIdentification.cpp`, which relies on acceleration measurements.

## Getting Started

Users should first take a look at the documentation of `TrajectoryData` class, where the trajectory data is imported and stored.

More documentation on python bindings will be coming!

## References
[1] C. Rucker and P. M. Wensing, "Smooth Parameterization of Rigid-Body Inertia," in IEEE Robotics and Automation Letters, vol. 7, no. 2, pp. 2771-2778, April 2022, doi: 10.1109/LRA.2022.3144517.


[2] K. Park and Y. Choi, "System identification method for robotic manipulator based on dynamic momentum regressor," 2016 12th IEEE International Conference on Control and Automation (ICCA), Kathmandu, Nepal, 2016, pp. 755-760, doi: 10.1109/ICCA.2016.7505369.