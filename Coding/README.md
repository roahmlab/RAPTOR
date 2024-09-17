# Coding Manual

This repo is coded in the KISS (Keep it stupid and simple) principle. 
There're private functions or members in any classes.
Any results are stored as class members so that they can directly accessed by others.

Readers can refer to `Examples/Kinova/include/KinovaOptimizer.h` and `Examples/Kinova/include/KinovaOptimizer.cpp` as a starting point to understand the structure.

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

Specifically, the core function of all trajectory classes is `compute(const Eigen::VectorXf& z, bool compute_derivatives)`.
This function takes in the decision variable `z` of the optimization problem and updates the class member `q`, `q_d`, `q_dd`.

 - `q` is an Eigen array of Eigen vectors.
 - `q(i)` stores all joint positions at time instance `i`.
 - `q(i)(j)` stores the position of joint `j` at time instance `i`.
 
 If `compute_derivatives` is true, then the function will also computes the gradient of `q`, `q_d`, `q_dd` with respect to the decision variable `z`, which are stored in `pq_pz`, `pq_d_pz`, and `pq_dd_pz`.

 - `pq_pz` is an Eigen array of Eigen matrices.
 - `pq_pz(i)` is a Eigen matrix of the jacobians of all joint positions with respect to the decision variable `z` at time instance `i`.

 These gradients are later taken in by other classes to eventually compute the gradients of the constraints or costs in the optimization.

For more information, please refer to [README](../Trajectories/README.md) in `Trajectories/` folder.

 ### Kinematics & Dynamics
For some of the constraints, such as end effector position, or torque limits, we need to compute the forward kinematics and inverse dynamics of the robot.
The core function is also called `compute(const Eigen::VectorXf& z, bool compute_derivatives)`

 - The forward kinematics or the inverse dynamics class usually requires a shared pointer of a trajectory class in the constructor.
 - In `compute`, the `compute` of the trajectory class will first be called, to update the trajectories of all joints.
 - Using `q`, `q_d`, `q_dd` in the trajectory class, we will compute the forward kinematics of each joints or the torque requried, depending on which kinematics/dynamics class it is.
 - If `compute_derivatives` is true, we will use `pq_pz`, `pq_d_pz`, `pq_dd_pz` in the trajectory class to compute the gradient of the forward kinematics or inverse dynamics with respect to the decision variable `z`.
 - The results will be stored inside the class for the constraints classes to access.

For more information, please refer to [README](../KinematicsDynamics/README.md) in `KinematicsDynamics/` folder. 

 ### Constraints
This class implements different types of constraints that will be later put into the optimization.
Specifically, we would like to provide the formulation `g_lb <= g(z) <= g_ub`.
The core function is called `compute(const Eigen::VectorXf& z, bool compute_derivatives)`. 
The constraint class usually requires a shared pointer of a trajectory class in the constructor, and a shared pointer of a forward kinematics and inverse dynamics class. 
The results are stored in `g`, which is a class member, for the optimizer class to access.
The gradient is stored in `pg_pz`.

Another core function is called `compute_bounds()`.
This is called by the optimizer class at the very beginning of the optimization.
This updates `g_lb` and `g_ub`, which are also class members.

There's also a function called `print_violation_information`, which will be called by the optimizer class at the end of the optimization, to print relevant information for tuning and debugging.

For more information, please refer to [README](../Constraints/README.md) in `Constraints/` folder. 

 ### Optimizer 
This class serves as an interface to an open-source nonlinear solver `Ipopt`.
We only cares about filling in the following functions:

 - `get_nlp_info`: initialize the number of the decision variable `z`, and the number of constraints.
 - `get_bounds_info`: initialize the bounds of the decision variable `z`, and the bounds of the constraints (where it will call `compute_bounds()` in each of the constraints class).
 - `eval_f`: provide the cost function value.
 - `eval_grad_f`: provide the gradient of the cost function.
 - `eval_g`: provide the constraints values.
 - `eval_jac_g`: provide the gradient (jacobian) of the constraints.
 - `eval_h`: provide the hessian of the entire Lagrange multiplier.
 - `finalize_solution`: `Ipopt` will return a solution here (may be infeasible or not reach the constraint violation tolerance yet), but you can check things here if you want to take this as a valid solution.


