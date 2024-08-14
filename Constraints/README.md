# Constraints

This class implements different types of constraints that will be later put into the optimization.
Specifically, we would like to provide the formulation `g_lb <= g(z) <= g_ub`.
The core function is called `compute(const Eigen::VectorXd& z, bool compute_derivatives)`. 
The constraint class usually requires a shared pointer of a trajectory class in the constructor, and a shared pointer of a forward kinematics and inverse dynamics class. 
The results are stored in `g`, which is a class member, for the optimizer class to access.
The gradient is stored in `pg_pz`.

Another core function is called `compute_bounds()`.
This is called by the optimizer class at the very beginning of the optimization.
This updates `g_lb` and `g_ub`, which are also class members.

There's also a function called `print_violation_information`, which will be called by the optimizer class at the end of the optimization, to print relevant information for tuning and debugging.

### JointLimits

Constrain robot joint angles by user-defined lower bounds and upper bounds.

### VelocityLimits

Constrain robot joint velocities by user-defined upper bounds.

### TorqueLimits

Constrain the torque of each joint within the limits.
We compute the torque using inverse dynamics implemented by [pinocchio](https://stack-of-tasks.github.io/pinocchio/).

### KinematicsConstraints

Constrain the position and the orientation of a specific joint at a specific time instance at a desired position and orientation.
This [website](https://wang-yimu.com/introduction-to-optimization-on-manifolds/) provides a pretty good introduction on how optimization should be done for this type of problem.
We have not implemented the axis-angle representation yet.
Instead, we constrain all 9 elements of the rotation matrix.
For axis-angle representation (SE(3)) of the constraints and the jacobian, please refer to this [paper](https://arxiv.org/abs/1606.05285).

### RectangleSurfaceContactConstraints

Assuming that an object is attached on the end effector of the robot.
The contact surface is flat and rectangle.
We need to satisfy certain constraints so that the object keeps a rigid contact with the robot.
We first need to access the contact wrench between the object and the robot.
The contact wrench contains 6 elements (`fx`, `fy`, `fz`, `mx`, `my`, `mz`).
`fx`, `fy`, `fz` represents the translational contact force in the local frame.
`mx`, `my`, `mz` represents the rotational contact moment in the local frame.
There are several conditions for the object to stay on the end effector:
