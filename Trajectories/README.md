# Trajectories

## Trajectories Class Overview

This folder contains various trajectory classes that take the decision variable `z` from the optimization problem and compute the corresponding trajectories for each joint `j` of the robot.

By trajectories, we refer to the:
- Positions, stored in public class member `q`, 
- Velocities, stored in public class member `q_d`, 
- Accelerations stored in public class member `q_dd`, 

over a series of predefined and fixed time instances.

The most common example is a uniform distribution of time instances over a fixed interval `[0, T]`, where `T` represents the total trajectory duration.

Another choice is Chebyshev nodes, to minimize the oscillations at the beginning and the end of trajectory.

### Fixed vs. Adaptive Time Instances

Currently, both trajectory duration (`T`) and the time distribution are fixed. However, in research literature, some approaches incorporate:
- Variable trajectory duration as an optimization decision variable to solve problems like minimum-time motion planning, while satisfying torque or velocity constraints.
- Adaptive time distributions, commonly explored in direct collocation methods and spectral methods for improved numerical performance.

For now, RAPTOR does not focus on these adaptive methods.

### Core Functions
Specifically, the core function of all trajectory classes is 

```C++
compute(const Eigen::VectorXd& z, bool compute_derivatives, bool compute_hessian)
```

This function takes the decision variable `z` from the optimization problem and updates the public class members `q`, `q_d`, and `q_dd`, which are later utilized by the `Constraints` class or the `Costs` class to compute the gradients of constraints and cost functions, facilitating efficient trajectory optimization.

#### Trajectory Representation
- `q` is an Eigen array of Eigen vectors, representing joint positions over time.
- `q(i)`: Stores all joint positions at time instance `i`.
- `q(i)(j)`: Stores the position of joint `j` at time instance `i`.

#### Gradient Computation
If `compute_derivatives` is `true`, the function also computes the gradients of `q`, `q_d`, and `q_dd` with respect to `z`. These are stored as:

- `pq_pz`: Eigen array of Eigen matrices, storing position gradients.
- `pq_pz(i)`: Jacobian matrix of all joint positions with respect to `z` at time instance `i`.

#### Hessian Computation
If `compute_hessian` is `true`, the function also computes the Hessian of `q`, `q_d`, and `q_dd` with respect to `z`. These are stored as:

- `pq_pz_pz`: 2D Eigen array of Eigen matrices, storing second-order derivatives.
- `pq_pz_pz(j, i)`: A square Eigen matrix representing the Hessian of the position of joint `j` at time instance `i` with respect to `z`.

## Trajectory Types

<!-- ### TrajectoryGroup class

`TrajectoryGroup` is a special type of the class that contains multiple trajectories.
For example, for multiple-step humanoid gait optimization, you can assign multiple trajectories in the group for each walking step of the humanoid robot and optimize all of these steps simultaneously.
For manipulation tasks, like the robotic arm picks one object and places it at another place, you can also assign one trajectory for pick motion and another trajectory for place motion, so that they can be optimized at the same time. -->

### Plain

This class implements a very naive "trajectory", which is just a point in the configuration space, in other words, the joint positions at one time instance.
The velocity and the acceleration are by default zero since there's no actual movement.

The number of decision variables (`varLength`): number of joints.

### TrajectoryData

This class does not define any trajectories, but load the time, joint positions and joint velocities from a file, which is used to load hardware data for system identification.
Otherwise, it generates random joint trajectories for users.
The file needs to contain a data matrix that
 - the first column is the time.
 - the second column to the n + 1 th column are the joint positions.
 - the n + 2 th column to the 2n + 1 th column are the joint velocities.
 - (optional) the 2n + 2 th column to the 3n + 1 th column are the joint accelerations or the joint torques.

Note that this class does not implement anything in `compute` function.
It does not fill in derivatives or hessians as well.

The number of decision variables (`varLength`): 0.

### Polynomials

This class implements a polynomial representation of the trajectory, where the coefficients of the polynomial are decision variables.

The number of decision variables (`varLength`): (number of joints) * (polynomial degree + 1).

### BezierCurves

This class implements a [Bezier curve](https://en.wikipedia.org/wiki/B%C3%A9zier_curve) representation of the trajectory, where the Bezier coefficients are decision variables.
Note that the Bezier curve is only defined on interval [0,1] originally.
We scale the curve to the interval [0,T] where T is the duration of the trajectory and a constant positive number.

You can use function `constrainInitialPosition(const VecX& q0_input)` to directly constrain the initial position of the trajectory to a specific vector.
In this instance, there will be one less decision variable for each joint.

Similarly, you can use function `constrainInitialVelocity(const VecX& q_d0_input)` to directly constrain the initial velocity of the trajectory to a specific vector.
In this instance, there will be one less decision variable for each joint.
In other words, if you use both `constrainInitialPosition` and `constrainInitialVelocity`, there will be two less decision variables.

The number of decision variables (`varLength`): 
 - (number of joints) * (Bezier curve degree + 1).
 - (number of joints) * (Bezier curve degree). (if either `constrainInitialPosition` or `constrainInitialVelocity` is activated)
 - (number of joints) * (Bezier curve degree - 1). (if both `constrainInitialPosition` and `constrainInitialVelocity` are activated)

### PiecewiseBezierCurves

This class implements a series of Bezier curves so that they are piecewise continuous on its second-order derivative.

Note that the Bezier curve has the following property:

 - The value at the beginning is equal to the first Bezier coefficient.
 - The first-order derivative at the beginning only depends on the first, and the second Bezier coefficient.
 - The second-order derivative at the beginning only depends on the first, the second, and the third Bezier coefficient.
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

In this class, you can also provide initial position (`q0_input` in the constructor inputs) so that the very beginning of the entire series of Bezier curves starts from this initial position with 0 velocity and 0 acceleration (since the robot usually starts from a static position).
Note that this initial position is not part of the decision variables and is treated as constant.
The same things hold for end position (`qT_input` in the constructor inputs), which is another optional input so that the robot stops at this position 0 velocity and 0 acceleration.
The motivation here is to directly constrain the start and the end of the entire trajectory using the property of Bezier curves, so that we have less decision variables and less constraints.

The number of decision variables (`varLength`): 
 - (number of joints) * degree * 3. (if both `q0_input` and `qT_input` are user-specified vectors in the constructor)
 - (number of joints) * (degree * 3 + 1). (if either `q0_input` or `qT_input` is a user-specified vector in the constructor)
 - (number of joints) * (degree * 3 + 2). 

### ArmourBezierCurves

This class is inherited from BezierCurves and describes a 5th-order Bezier curve that starts from specific position, velocity, acceleration, and ends at a desired position (which is the only decision variable here) with 0 velocity and 0 acceleration.
That means there's only 1 decision variable for each joint of the robot.
It is used in our previous work [ARMOUR](https://roahmlab.github.io/armour/).

The number of decision variables (`varLength`): number of joints

### FourierCurves

This class represents trajectories using a Fourier series expansion. Different from the previous work, **RAPTOR** expresses joint accelerations as a sum of sinusoidal basis functions.

Consider `N` as the number of harmonics used, also known as 'degree' here, the joint acceleration is represented as:

```math
\ddot{q}(t) = a_0 + \sum_{k=1}^{N} \left( a_k \cos(k\omega t) + b_k \sin(k\omega t) \right)
```
where:
- `a_0` is the constant offset,
- `a_k` and `b_k` are Fourier coefficients,
- `\omega` is the base/fundamental frequency.

All of the variables above are decision variables for this trajectory class.

Similar to the `BezierCurves` class, it provides options to optimize the initial position or the initial velocity or not.
Users can provide `q0_input` in the constructor as an input to start the trajectory at a specific position, otherwise, **RAPTOR** will treat the initial position also as decision variables.
This is similar to the initial velocity as well (`q_d0_input` in the constructor inputs).

The number of decision variables (`varLength`): 
 - (number of joints) * (2 * degree + 2). (if both `q0_input` and `qT_input` are user-specified vectors in the constructor)
 - (number of joints) * (2 * degree + 3). (if either `constrainInitialPosition` or `constrainInitialVelocity` is activated)
 - (number of joints) * (2 * degree + 4). 

Note that the Hessian of this trajectory class is not implemented.

### FixedFrequencyFourierCurves

Similar to the `FourierCurves` class, this class also implements a Fourier series expansion but with a fixed base/fundamental frequency.
In other words, you have to provide `\omega` as an input.

The number of decision variables (`varLength`): 
 - (number of joints) * (2 * degree + 1). (if both `q0_input` and `qT_input` are user-specified vectors in the constructor)
 - (number of joints) * (2 * degree + 2). (if either `constrainInitialPosition` or `constrainInitialVelocity` is activated)
 - (number of joints) * (2 * degree + 3).