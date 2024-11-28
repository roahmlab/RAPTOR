# Trajectories

## Overview

### Trajectories class

This folder contains a variety of trajectory classes that take in the decision variable `z` of the optimization problem and returns the trajectories of each joint `j` of the robot.
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

If `compute_hessian` is true, then the function will also computes the hessian of `q`, `q_d`, `q_dd` with respect to the decision variable `z`, which are stored in `pq_pz_pz`, `pq_d_pz_pz`, and `pq_dd_pz_pz`.
 - `pq_pz_pz` is a 2D Eigen array of Eigen matrices.
 - `pq_pz_pz(j, i)` is a square Eigen matrix of the hessian of the position with respect to the decision variable `z` for joint `j` at time instance `i`. 

These gradients or hessians are later taken in by other classes to eventually compute the gradients of the constraints or costs in the optimization.

### TrajectoryGroup class

`TrajectoryGroup` is a special type of the class that contains multiple trajectories.
For example, for multiple-step humanoid gait optimization, you can assign multiple trajectories in the group for each walking step of the humanoid robot and optimize all of these steps simultaneously.
For manipulation tasks, like the robotic arm picks one object and places it at another place, you can also assign one trajectory for pick motion and another trajectory for place motion, so that they can be optimized at the same time.

## Introduction to Each Trajectories

### Plain

This class implements a very naive "trajectory", which is just the joint positions at one time instance.
The velocity and the acceleration are by default zero since there's no actual movement.

### TrajectoryData

This class does not define any trajectories, but load the time, joint positions and joint velocties from a file.
Otherwise, it generates random joint trajectories for users.
The file needs to contain a data matrix that
 - the first column is the time.
 - the second column to the n + 1 th column are the joint positions.
 - the n + 2 th column to the 2n + 1 th column are the joint velocities.

Note that this class does not implement anything in `compute` function.
It does not fill in derivatives or hessians as well.

### Polynomials

This class implements a polynomial representation of the trajectory, where the coefficients of the polynomial are decision variables.

### BezierCurves

This class implements a [Bezier curve](https://en.wikipedia.org/wiki/B%C3%A9zier_curve) representation of the trajectory, where the Bezier coefficients are decision variables.
Note that the Bezier curve is only defined on interval [0,1] originally.
We scale the curve to the interval [0,T] where T is the duration of the trajectory and a constant positive number.

### PiecewiseBezierCurves

This class implements a series of Bezier curves so that they are piecewise continuous on its second-order derivative.

Note that the Bezier curve has the following property:

 - The value at the beginning is equal to the first Bezier coefficient.
 - The first-order derivative at the beginning only depends on the first, and the second Bezier coefficient.
 - The second-order derivative at the beginning only depends on the first, the second, and the thrid Bezier coefficient.
 - The same property hold for the value, the first-order derivative, and the second-order derivative at the end, which only depends on the last first, the last second, and the last third Bezier coefficient.

As a result, here we define a 5th-order Bezier curve which has 6 Bezier coefficients.
We consider the position, the velocity, and the acceleration at the beginning and at the end, which are also 6 numbers that uniquely determine the 6 Bezier coefficients.
We simplify the notation here and define the position, the velocity, and the acceleration at the beginning as `y_{i}` and at the end as `y_{i+1}`.
`y_{i}` and `y_{i+1}` are considered as decision variables.
In this class, we will transform `y_{i}` and `y_{i+1}` to corresponding Bezier coefficients so that we compute the trajectory at arbitrary time instances.

Now let's consider a series of Bezier curves.
Let's say the `i` th Bezier curve has `y_{i}` as the position, velocity, acceleration at the beginning and `y_{i+1}` as the position, velocity, acceleration at the end.
The `i+1` th Bezier curve has `y_{i+1}` as the position, velocity, acceleration at the beginning and `y_{i+2}` as the position, velocity, acceleration at the end.
The `i` th Bezier curve and the `i+1` th Bezier curve thus are 2nd-order continuously differentiable (continuous on its 2nd order derivative) at their connection point, since they share the same positon, velocity, and acceleration here.

As a result, we consider a series of `y_{i}` with `i=1,...N`.
We then can define `N-1` Bezier curves that are connected with each other in series.

In this class, you can also provide initial position so that the very beginning of the entire series of Bezier curves starts from this initial position with 0 velocity and 0 acceleration (since the robot usually starts from a static position).
Note that this initial position is not part of the decision variables and is treated as constant.
The same things hold for end position, which is another optional input so that the robot stops at this position 0 velocity and 0 acceleration.
The motivation here is to directly constrain the start and the end of the entire trajectory using the property of Bezier curves, so that we have less decision variables and less constraints.

### ArmourBezierCurves

This class is inherited from BezierCurves and describes a 5th-order Bezier curve that starts from specific position, velocity, acceleration, and ends at a desired position (which is the only decision variable here) with 0 velocity and 0 acceleration.
That means there's only 1 decision variable for each joint of the robot.
It is used in our previous work [ARMOUR](https://roahmlab.github.io/armour/).