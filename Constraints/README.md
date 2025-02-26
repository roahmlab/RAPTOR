# Constraints

## Class Overview
This class implements various types of constraints that are later incorporated into the optimization process. The constraints are formulated as:
```C++
g_lb <= g(z) <= g_ub
```

The primary function is 
```C++
compute(const Eigen::VectorXd& z, bool compute_derivatives, bool compute_hessian)
``` 
This function evaluates the constraint function `g(z)`.
 - The computed constraint values are stored in an Eigen vector `g`, a public class member, allowing the optimizer to access them from outside.
 - If compute_derivatives is true, the gradient (the jacobian matrix) is computed and stored in an Eigen matrix `pg_pz`.
 - If compute_hessian is true, the Hessian is computed and stored in `pg_pz_pz`. Note that `pg_pz_pz` is an 1-D Eigen array of Eigen matrices, where `pg_pz_pz(i)` stores the Hessian matrix of the `i`th constraint.

Another core function is called 
```C++
compute_bounds()
```
This is called by the optimizer class at the very beginning of the optimization.
This updates `g_lb` and `g_ub`, which are also public class members, allowing the optimizer to access them from outside.

There's also a function called 
```C++
print_violation_information()
```
which will be called by the optimizer class at the end of the optimization, to print relevant information for tuning and debugging.

## Constraint Types

### JointLimits

Constrain robot joint angles by user-defined lower bounds and upper bounds.

Hessian implementation: **Yes**

### ConstrainedJointLimits

A special version of `JointLimits` for constrained systems (systems with kinematics constraints, such as [Digit-v3](../Examples/Digit/), which incorporates closed-loop kinematics constraints at its knees and ankles)

Hessian implementation: **No**

### VelocityLimits

Constrain robot joint velocities by user-defined upper bounds.

Hessian implementation: **Yes**

### AccelerationLimits

Constrain robot joint accelerations by user-defined upper bounds.

Hessian implementation: **Yes**

### TrajectoryTerminalConstraints

Constrain robot position, velocity, or acceleration at the end of the trajectory by user-defined values.

### TorqueLimits

Constrain the torque of each joint within the limits.
We compute the torque using inverse dynamics implemented by [pinocchio](https://stack-of-tasks.github.io/pinocchio/).

Hessian implementation: **Yes**

### KinematicsConstraints

Constrain the position and the orientation of a specific joint at a specific time instance at a desired position and orientation.
This [website](https://wang-yimu.com/introduction-to-optimization-on-manifolds/) provides a pretty good introduction on how optimization should be done for this type of problem.
We have not implemented the axis-angle representation yet.
Instead, we constrain all 9 elements of the rotation matrix.
For axis-angle representation (SE(3)) of the constraints and the jacobian, please refer to this [paper](https://arxiv.org/abs/1606.05285).

Hessian implementation: **Yes**

### RectangleSurfaceContactConstraints

Enforcing contacts between the robot and one object or the ground, where the contact surface is **rectangle**.

For example, the object to manipulate should be attached on the end effector of the robot, or the robot stays on the ground on its stance foot.
We need to satisfy certain constraints so that the object keeps a rigid contact with the robot.
We first need to access the contact wrench between the object and the robot.
The contact wrench contains 6 elements (`fx`, `fy`, `fz`, `mx`, `my`, `mz`).
`fx`, `fy`, `fz` represents the translational contact force in the local frame.
`mx`, `my`, `mz` represents the rotational contact moment in the local frame.
There are several conditions for the object to stay on the end effector:
```math
\begin{align}
    0 &\leq \lambda_{st}^{fz}(t) \\
    \sqrt{(\lambda_{st}^{fx}(t))^2 + (\lambda_{st}^{fy}(t))^2} &\leq \mu \lambda_{st}^{fz}(t) \\
    \lambda_{st}^{mz}(t) &\leq \gamma \lambda_{st}^{fz}(t) \\
    -\frac{1}{2}l_a\lambda_{st}^{fz}(t) \leq \lambda_{st}^{mx}(t) &\leq \frac{1}{2}l_a\lambda_{st}^{fz}(t) \\
    -\frac{1}{2}l_b\lambda_{st}^{fz}(t) \leq \lambda_{st}^{my}(t) &\leq \frac{1}{2}l_b\lambda_{st}^{fz}(t)
\end{align}
```

- (1) states that the contact force should always be non-negative.
- (2) states that the translational friction force should stay inside the friction cone.
- (3) states that the rotational friction force should stay inside the friction cone.
- (4) and (5) state the ZMP constraints, that the object should not roll over one of the edge of the rectangle contact surface.

Hessian implementation: **No**

### CircleSurfaceContactConstraints

Enforcing contacts between the robot and one object or the ground, where the contact surface is **circle**.
The constraints are formulated similarly to `RectangleSurfaceContactConstraints`, but with a much simpler ZMP constraints formulation.

Hessian implementation: **No**

## A Special Type of "Constraint": CollisionAvoidance

`CollisionAvoidance` is actually a separate class, which is not inherited from `Constraints` class.
It aims to compute the minimum distance between the environment and a point.
Given the gradient and the Hessian of the point (with respect to certain decision variables), it also computes the gradient and the Hessian of the distance.
It can be used in constraint classes that enforce collision avoidance, in other words, the distance between the environment and a specific point on the robot greater than a threshold.
An example of such constraint class can be found in `KinovaCustomizedConstraints`.

In **RAPTOR**, we have only implemented a class called `BoxCollisionAvoidance` inherited from `CollisionAvoidance`.
The `BoxCollisionAvoidance` class takes a vector of boxes representing objects in the environment and computes the minimum distance between a specified point and all the boxes.
