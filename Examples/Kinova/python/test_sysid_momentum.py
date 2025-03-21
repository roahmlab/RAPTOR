import numpy as np
import pinocchio as pin
import scipy.io as sio
import sys

from kinova_dynamics import integrate

sys.path.append("/workspaces/RAPTOR/build/lib")
import end_effector_sysid_momentum_nanobind

def desired_trajectory(t):
    """
    Computes the desired trajectory for a 7-DOF system based on sinusoidal functions.

    Parameters:
    t (float): The time variable.

    Returns:
    tuple: A tuple containing:
        - qd (numpy.ndarray): The desired joint positions, a 7-element array where each element is sin(t).
        - qd_d (numpy.ndarray): The desired joint velocities, a 7-element array where each element is cos(t).
        - qd_dd (numpy.ndarray): The desired joint accelerations, a 7-element array where each element is -sin(t).
        
    Note:
        This trajectory is not an exciting one, so the results may not be good on hardware
    """
    qd = np.sin(t) * np.ones(7)
    qd_d = np.cos(t) * np.ones(7)
    qd_dd = -np.sin(t) * np.ones(7)
    return qd, qd_d, qd_dd

def controller(q, v, qd, qd_d, qd_dd):
    """
    Computes the control torque for a system using a PD control law.

    Parameters:
        q (float or array-like): Current position of the system.
        v (float or array-like): Current velocity of the system.
        qd (float or array-like): Desired position of the system.
        qd_d (float or array-like): Desired velocity of the system.
        qd_dd (float or array-like): Desired acceleration of the system (not used in this implementation).

    Returns:
        float or array-like: Control torque to be applied to the system.
    """
    
    kp = 1000
    kd = 100
    tau = kp * (qd - q) + kd * (qd_d - v)
    return tau

def main():
    # initialization for simulation and data collection
    urdf_filename = "/workspaces/RAPTOR/Robots/kinova-gen3/gen3_2f85_fixed.urdf"
    model = pin.buildModelFromUrdf(urdf_filename)
    
    q0 = np.zeros(model.nq)
    v0 = np.zeros(model.nv)
    
    dt = 1e-4 # 0.1 ms data measurement loop
    ts_sim = np.arange(0, 5, dt) # 5 seconds simulation
    
    # simulate the robot dynamics using ode solver
    # track the desired trajectory using the controller
    qs, vs, taus = integrate(model, ts_sim, np.concatenate([q0, v0]), desired_trajectory, controller)
    
    # save the simulation results as a text file
    traj_data = np.concatenate([ts_sim[:,None], qs, vs, taus], axis=1)
    np.savetxt("traj_data.csv", traj_data, delimiter=" ")
    
    # initialization for system identification
    H = 10 # forward integration horizon
    
    # motor friction parameters
    # [ --- static friction                 ---
    #   --- viscous friction (damping)      ---
    #   --- transmission inertia (armature) ---
    #   --- offset                          ---]
    # in this particular example, motor friction parameters are disabled (set to zero)
    friction_parameters = np.zeros(28)

    sysid_solver = end_effector_sysid_momentum_nanobind.EndEffectorIdentificationMomentumPybindWrapper(
        urdf_filename, 
        friction_parameters, 
        H,                   # forward integration horizon
        'Second',            # format of time data, 'Second' or 'Nanosecond'
        True                 # display information or not
    )
    
    sysid_solver.add_trajectory_file("traj_data.csv")

    sysid_solver.set_ipopt_parameters(
        1e-14,          # tol
        5.0,            # max_wall_time
        5,              # print_level
        500,            # max_iter
        "adaptive",     # mu_strategy
        "ma86",         # linear_solver
        False           # gradient_check
    )

    theta_solution, theta_uncertainty = sysid_solver.optimize()
    
    print("end effector inertial parameters:")
    print("solution:\n", theta_solution)
    print("uncertainty:\n", theta_uncertainty)
    print("groundtruth:\n", model.inertias[-1].toDynamicParameters())

if __name__ == "__main__":
    main()
