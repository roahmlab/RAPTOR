# Kinova-gen3 Examples

This folder contains multiple optimization examples related to the Kinova-gen3, a 7-degree-of-freedom robotic manipulator. A Python interface is also provided to enable the use of these examples in Python.

## CollisionAvoidanceTrajectory

This folder contains two examples of trajectory generation that satisfy the following constraints over a predefined set of time instances:

- Joint limits
- Velocity limits
- Torque limits
- Collision avoidance with a series of box obstacles

In addition to ensuring safety through these constraints, the following optimization formulations are considered:

### KinovaOptimizer

The trajectory is represented using `ArmourBezierCurves`, a degree-5 Bézier curve with zero terminal velocity and acceleration. The cost function is defined as the distance between the midpoint of the trajectory and a user-defined goal configuration.

This formulation can be seen as a discrete version of [ARMOUR](https://roahmlab.github.io/armour/), a receding-horizon trajectory optimization framework. Given a sequence of waypoints (e.g., from a high-level path planner like RRT), ARMOUR optimizes a degree-5 Bézier curve at each iteration to move the robot closer to the next waypoint while respecting all constraints that the high-level planner might not account for.

### KinovaLongerHorizonOptimizer

This trajectory optimizer finds a safe trajectory that connects a user-defined start configuration to a goal configuration. The trajectory is represented using `PiecewiseBezierCurves`.

The cost function is customizable, allowing optimization for various objectives, such as:
- Minimizing overall torque consumption
- Minimizing path length
- Minimizing jerk for smoother and more energy-efficient motions

## CollisionAvoidanceInverseKinematics

This folder contains an example of solving inverse kinematics with collision avoidance constraints for a Kinova-gen3 equipped with a gripper. The following constraints are enforced:

- Joint limits.
- End-effector kinematic constraints, including position and orientation constraints.
- Collision avoidance with a series of box obstacles.

Since inverse kinematics deals with finding a single joint configuration rather than a full trajectory, the configuration is represented using `Plain`.

The cost function minimizes the distance from the initial guess. 
Conceptually, this can be interpreted as:

- The robot is at a given configuration.
- It must move to a new configuration that satisfies specific position and orientation constraints while avoiding obstacles.
- The cost function ensures that the final configuration remains as close as possible to the initial guess.

A key limitation of this formulation is its sensitivity to the initial guess. If the initial guess is poor, the optimization may struggle to find a feasible solution. However, this issue is less frequent for the Kinova-gen3, as it has 7 degrees of freedom (DOF), providing more flexibility in satisfying constraints.
For 6 DOF robots, users may need to consider alternative formulations to improve robustness.

## Armour

This folder contains an implementation that integrates both [ARMOUR](https://roahmlab.github.io/armour/) and [WAITR](https://roahmlab.github.io/waitr-dev/).

Building on the discussion in CollisionAvoidanceTrajectory, traditional trajectory optimization formulations typically enforce constraints only at discrete time instances. 
This means that safety between time steps is not guaranteed.

Our previous work, [ARMTD](https://arxiv.org/abs/2002.01591), first addressed this issue for robotic manipulators in the context of collision avoidance. 
This was further extended in [ARMOUR](https://roahmlab.github.io/armour/), which not only ensured continuous-time safety in collision avoidance but also enforced torque limits.
[WAITR](https://roahmlab.github.io/waitr-dev/) builds upon these concepts by providing safety guarantees for contact constraints. 
Specifically, it ensures stability when the robot gripper holds a plate with an object not rigidly attached, a problem commonly known as the waiter motion problem.

This folder provides an implementation that integrates all of our previous work. 
By modifying the configuration YAML file, users can solve different types of problems based on their specific requirements.
More details can be found in the [Armour](../Kinova/Armour/) folder.

The Hessian of this optimization problem is not provided.

## SystemIdentification

This folder contains several examples for system identification of Kinova-gen3.
There are two main topics in the context of system identification:

### Exciting Trajectories

Exciting trajectories are specifically designed motion paths that maximize the information gained about the system's dynamics. In the context of system identification, these trajectories are crucial because they help in accurately estimating the parameters of the Kinova-gen3 robotic manipulator.

By carefully designing these trajectories, we can ensure that the collected data is rich in information, which leads to more precise and reliable system models. These models are essential for various applications, including control design, simulation, and optimization.

The examples in the [ExcitingTrajectories](../Kinova/SystemIdentification/ExcitingTrajectories/) folder demonstrate how to generate and utilize exciting trajectories to improve the accuracy of system identification for the Kinova-gen3 while ensuring safety, such as joint limits, torque limits, or collision avoidance.
More details can be found in the folder.

### Parameter Identification

The [ParametersIdentification](../Kinova/SystemIdentification/ParametersIdentification/) folder contains the following examples:

#### Friction Identification
In the context of robotic manipulators, motor friction can significantly affect the accuracy and performance of the system. 
To model the friction in the motors of the Kinova-gen3, we can use a simpler model that includes offset, static friction, damping, and armature (transmission inertia) components. 
The mathematical model for the motor friction `\tau_f` can be expressed as:

```math
\tau_f(\dot{q}, \ddot{q}) = \beta + F_c \cdot \text{sign}(\dot{q}) + F_v \cdot \dot{q} + I_a \cdot \ddot{q}
```
where:
- `\beta` is the offset torque, which accounts for any constant bias in the system.
- `F_c` is the static friction.
- `F_v` is the damping coefficient, which represents the torque proportional to the joint velocity `\dot{q}`.
- `I_a` is the armature (transmission inertia), which represents the torque proportional to the joint acceleration `\ddot{q}`.

By assuming that the inertial parameters of the robot have been provided by the manufacturer (the URDF of the robot), we can identify these motor friction parameters by following certain exciting trajectories and recording the joint positions, velocities and applied torques.

#### End Effector Identification

Assuming that the robot's inertial parameters and motor friction parameters are known, accurate dynamic modeling can be achieved for the robot itself. However, in manipulation tasks, when the robot picks up an unknown object, its overall dynamics change.
To ensure precise control and stability, it becomes crucial for the robot to identify the inertial parameters of the object being held. 
More specifically, this involves estimating the combined inertial properties of the end effector and the object as a single entity. 
This estimation allows for more accurate force and motion planning, improving performance in tasks such as object transport, placement, and interaction with the environment.

The [ParametersIdentification](../Kinova/SystemIdentification/ParametersIdentification/) folder contains two different methods to estimate the end effector inertial parameters.
One is based on the inverse dynamics regressor which requires estimation of joint acceleration.
The other is based on the system momentum regressor which does not require joint acceleration.





