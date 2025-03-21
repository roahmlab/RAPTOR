import numpy as np
import pinocchio as pin
import scipy.io as sio
from scipy.signal import butter, filtfilt
import sys

from kinova_dynamics import integrate

sys.path.append("/workspaces/RAPTOR/build/lib")
import end_effector_sysid_nanobind

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

def central_difference_4th_order(t, velocity):
    """
    Apply a 4th order central difference method to estimate acceleration 
    from velocity data.

    Parameters:
        t (numpy.ndarray): A 1D array of time values. Assumes uniform time steps.
        velocity (numpy.ndarray): A 2D array of velocity values with shape (n, m),
                                  where n is the number of time steps and m is the 
                                  number of velocity components.

    Returns:
        tuple:
            - numpy.ndarray: A 2D array of estimated accelerations with shape 
                             (n-4, m), corresponding to the interior points 
                             where the 4th order central difference is applied.
            - float: The average time step size (dt) computed from the time array.

    Notes:
        - The function assumes uniform time steps in the input time array.
        - Boundary points (first two and last two rows) are excluded from the 
          output as they do not have enough neighbors for a 4th order difference.
        - Boundary accelerations are padded with NaN values internally but are 
          not included in the returned acceleration array.
    """
    dt = np.diff(t).mean()  # Assuming uniform time steps
    n = velocity.shape[0]
    acceleration = np.zeros((n, velocity.shape[1]))

    # Apply 4th order central difference formula for interior points
    for i in range(2, n - 2):
        acceleration[i, :] = (
            -velocity[i + 2, :]
            + 8 * velocity[i + 1, :]
            - 8 * velocity[i - 1, :]
            + velocity[i - 2, :]
        ) / (12 * dt)

    # For boundary points, pad with NaN/0
    acceleration[:2] = np.nan  # Not enough points for 4th order
    acceleration[-2:] = np.nan

    return acceleration[2:-2], dt

def butterworth_lowpass_filter(data, cutoff, fs, order=4):
    """
    Apply a Butterworth low-pass filter to the input data.
    This function applies a digital Butterworth low-pass filter to the given 
    multi-dimensional data. The filter is applied independently to each column 
    of the input data.
    Parameters:
        data (numpy.ndarray): A 2D array where each column represents a signal 
            to be filtered.
        cutoff (float): The cutoff frequency of the low-pass filter in Hz.
        fs (float): The sampling frequency of the input data in Hz.
        order (int, optional): The order of the Butterworth filter. Default is 4.
    Returns:
        numpy.ndarray: The filtered data with the same shape as the input data.
    Notes:
        - The function uses zero-phase filtering via `scipy.signal.filtfilt` 
          to avoid phase distortion.
        - Ensure that the cutoff frequency is less than half the sampling 
          frequency (Nyquist frequency) to avoid aliasing.
    """
    b, a = butter(order, cutoff, btype='low', analog=False, fs = fs)
    data_filtered = data
    
    for i in range(data.shape[1]):
        data_filtered[:, i] = filtfilt(b, a, data[:, i])
    return data_filtered

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
    
    # estimate acceleration using central difference method on velocity data
    accs, dt = central_difference_4th_order(ts_sim, vs)
    fs = 1 / dt
    
    # # filter acceleration data using a Butterworth low-pass filter
    # cutoff = 20 # Hz
    # accs_filtered = butterworth_lowpass_filter(accs, cutoff, fs)
    accs_filtered = accs # filtering is not needed in simulation, but really important for hardware data
    
    # save the simulation results as a text file
    traj_data = np.concatenate([ts_sim[:,None], qs, vs, taus], axis=1)
    traj_data_clipped = traj_data[2:-2, :]
    np.savetxt("traj_data.csv", traj_data_clipped, delimiter=" ") # time, position, velocity, torque go together
    np.savetxt("acceleration_filtered.csv", accs_filtered, delimiter=" ") # acceleration goes separately
    
    # initialization for system identification
    
    # motor friction parameters
    # [ --- static friction                 ---
    #   --- viscous friction (damping)      ---
    #   --- transmission inertia (armature) ---
    #   --- offset                          ---]
    # in this particular example, motor friction parameters are disabled (set to zero)
    friction_parameters = np.zeros(28)

    sysid_solver = end_effector_sysid_nanobind.EndEffectorIdentificationPybindWrapper(
        urdf_filename, 
        friction_parameters,
        'Second',            # format of time data, 'Second' or 'Nanosecond'
        True                 # display information or not
    )
    
    sysid_solver.add_trajectory_file(
        "traj_data.csv",
        "acceleration_filtered.csv")

    sysid_solver.set_ipopt_parameters(
        1e-14,          # tol
        5.0,            # max_wall_time
        5,              # print_level
        500,            # max_iter
        "adaptive",     # mu_strategy
        "ma86",         # linear_solver
        False           # gradient_check
    )

    theta_solution = sysid_solver.optimize()
    
    print("end effector inertial parameters:")
    print("solution:\n", theta_solution)
    print("groundtruth:\n", model.inertias[-1].toDynamicParameters())

if __name__ == "__main__":
    main()
