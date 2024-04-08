# Coding Manual

This repo is coded in the KISS (Keep it stupid and simple) principle. 
There're private functions or members in any classes.
Any results are stored as class members so that they can directly accessed by others.
There are several stages of implementation:

### Trajectories
This class takes in the decision variable `z` of the optimization problem and returns the trajectories of each joint `j` of the robot.
By trajectories, we mean the position, velocity, and the acceleration of each joint `j` over a series of time instances.
This series of time instances is predefined and fixed during the optimization process.
The most common example would be a uniform distribution over a fixed time interval `[0,T]`, where `T` is the duration of the trajectory.

Right now, this duration and the distribution of time instances are fixed.
There are papers that incorporate the duration as one of the optimization decision variable, in the favor of solving problems, for example, reaching some target with the minimum amount of time while satisfying the torque limits or velocity limits.
Adaptive changes over the distribution of time instances are also widely discussed in the context of direct collocation methods or spectral methods.
But these are not the focus of our method for now.

Specifically, the core function of all trajectory classes is `compute(const Eigen::VectorXd& z, bool compute_derivatives)`.
This function takes in the decision variable `z` of the optimization problem and updates the class member `q`, `q_d`, `q_dd`.

 - `q` is an Eigen array of Eigen vectors.
 - `q(i)` stores all joint positions at time instance `i`.
 - `q(i)(j)` stores the position of joint `j` at time instance `i`.
 
 If `compute_derivatives` is true, then the function will also computes the gradient of `q`, `q_d`, `q_dd` with respect to the decision variable `z`, which are stored in `pq_pz`, `pq_d_pz`, and `pq_dd_pz`.

 - `pq_pz` is an Eigen array of Eigen matrices.
 - `pq_pz(i)` is a Eigen matrix of the jacobians of all joint positions with respect to the decision variable `z` at time instance `i`.

 These gradients are later taken in by other classes to eventually compute the gradients of the constraints or costs in the optimization.

 ### Kinematics & Dynamics
For some of the constraints, such as end effector position, or torque limits, we need to compute the forward kinematics and inverse dynamics of the robot.
The core function is also called `compute(const Eigen::VectorXd& z, bool compute_derivatives)`

 - The forward kinematics or the inverse dynamics class usually requires a shared pointer of a trajectory class in the constructor.
 - In `compute`, the `compute` of the trajectory class will first be called, to update the trajectories of all joints.
 - Using `q`, `q_d`, `q_dd` in the trajectory class, we will compute the forward kinematics of each joints or the torque requried, depending on which kinematics/dynamics class it is.
 - If `compute_derivatives` is true, we will use `pq_pz`, `pq_d_pz`, `pq_dd_pz` in the trajectory class to compute the gradient of the forward kinematics or inverse dynamics with respect to the decision variable `z`.

 ### Constraints


 ### Optimizer 

