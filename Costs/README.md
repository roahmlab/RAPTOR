# Costs

## Class Overview
This class implements different types of cost functions that will be later put into the optimization.

The primary function is
```C++
compute(const Eigen::VectorXd& z, bool compute_derivatives, bool compute_hessian)
``` 

The objective function value is stored in a `double` variable `f`, which is a public class member, for the optimizer class to access from outside.
The gradient is stored in an Eigen vector `grad_f` and the hessian is stored in an Eigen matrix `hess_f`.

The cost class usually requires a shared pointer of a trajectory class in the constructor, or a shared pointer of a forward kinematics and inverse dynamics class depending on the purpose of the cost functions. 

## Cost Types

### MinimizeInitialVelocity
Minimize (the norm of) the initial velocity of the trajectory.

### MinimizeInitialAcceleration
Minimize (the norm of) the initial acceleration of the trajectory.

### MinimizePathLength
Minimize the length of the trajectory.
The trajectory length is computed by summing the norms of joint position differences between consecutive time steps:

```math
f = \sum_{t=1}^{N-1} \| q(t_{i+1}) - q(t_{i}) \|
```

### Minimize Jerk

This cost function minimizes the jerk, which is the third derivative of position with respect to time. Jerk minimization ensures smooth motion by reducing sudden changes in acceleration.
The cost function here is computed as:

```math
f = \sum_{i=1}^{N} \| \dddot{q}(t_i) \| \propto \sum_{t=1}^{N-1} \| \ddot{q}(t_{i+1}) - \ddot{q}(t_{i}) \|
```

### MinimizeTorque
This cost function minimizes the torque applied at each joint throughout the trajectory, promoting energy efficiency and smoother control, which is computed as:

```math
f = \sum_{i=1}^{N} \| \tau(t_i) \|
```

### MinimizePower
This cost function minimizes the power consumption during motion, encouraging energy-efficient actuation, which is computed as:

```math
f = \sum_{i=1}^{N} \| \tau(t_i) \dot{q}(t_i) \|
```

### EndEffectorRegressorConditionNumber

This cost function minimizes the condition number of the regressor matrix associated with the end effector inertial parameters. 
A well-conditioned regressor improves numerical stability and parameter estimation accuracy in system identification.

Here, we minimize the difference between the logarithms of the maximum and minimum singular values of the regressor matrix, which is equivalent to the condition number (the ratio between the maximum and minimum singular values of the regressor matrix).

Note that the Hessian of this cost is not implemented.