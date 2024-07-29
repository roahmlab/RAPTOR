import numpy as np
import build.oracle_nanobind as oracle

urdf_filename = "Robots/kinova-gen3/kinova.urdf"

planner = oracle.KinovaPybindWrapper(urdf_filename)

# obstacle information (xyz, rpy, size)
obstacles = np.array([[0, 0, 2, 0, 0, 0, 0.1, 0.1, 0.1],
                      [2, 0, 0, 0, 0, 0, 0.1, 0.1, 0.1],
                      [0, 2, 0, 0, 0, 0, 0.1, 0.1, 0.1]])

# trajectory information
q0 = np.array([0, 0.5, 0, -0.7, 0.5, 0, 0.5])
qd0 = np.array([0, -1, 0, 0, 0, -1, 0])
qdd0 = np.array([1, 0, 0, 1, 0, 0, 0])

duration = 2.0

# target information
q_des = np.array([1, 1, 1, 1, 1, 1, 1])
t_plan = 0.5 * duration

# ipopt settings
ipopt_tol = 1e-3
ipopt_obj_scaling_factor = 1.0
ipopt_max_wall_time = 0.5
ipopt_print_level = 5
ipopt_mu_strategy = "adaptive"
ipopt_linear_solver = "ma57"
ipopt_gradient_check = False

# buffer information
joint_limits_buffer = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
velocity_limits_buffer = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
torque_limits_buffer = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

planner.set_obstacles(obstacles)
planner.set_ipopt_parameters(ipopt_tol, \
                             ipopt_obj_scaling_factor, \
                             ipopt_max_wall_time, \
                             ipopt_print_level, \
                             ipopt_mu_strategy, \
                             ipopt_linear_solver,
                             ipopt_gradient_check)
planner.set_trajectory_parameters(q0, \
                                  qd0, \
                                  qdd0, \
                                  duration)
planner.set_buffer(joint_limits_buffer, \
                   velocity_limits_buffer, \
                   torque_limits_buffer)
planner.set_target(q_des, t_plan)

[k, ifFeasible] = planner.optimize()

print(k)