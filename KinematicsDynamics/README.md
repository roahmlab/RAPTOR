# Kinematics & Dynamics

For some of the constraints, such as end effector position, or torque limits, we need to compute the forward kinematics and inverse dynamics of the robot.
The core function is also called `compute(const Eigen::VectorXd& z, bool compute_derivatives, bool compute_hessian)`

 - The forward kinematics or the inverse dynamics class usually requires a shared pointer of a trajectory class in the constructor.
 - In `compute`, the `compute` of the trajectory class will first be called, to update the trajectories of all joints.
 - Using `q`, `q_d`, `q_dd` in the trajectory class, we will compute the forward kinematics of each joints or the torque requried, depending on which kinematics/dynamics class it is.
 - If `compute_derivatives` is true, we will use `pq_pz`, `pq_d_pz`, `pq_dd_pz` in the trajectory class to compute the gradient of the forward kinematics or inverse dynamics with respect to the decision variable `z`.
 - The results will be stored inside the class for the kinematics/dynamics classes to access.