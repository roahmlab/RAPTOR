import numpy as np
import pinocchio as pin
import pybullet as p
import math
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

def goal_distance(diff):
    """
    Computes the Euclidean distance of a vector after wrapping specific angles 
    to the range [-π, π].

    The function adjusts the angles at indices 0, 2, 4, and 6 of the input 
    vector `diff` to ensure they are within the range [-π, π]. It then calculates 
    the Euclidean norm of the resulting vector.

    Parameters:
        diff (numpy.ndarray): A vector containing angular differences.

    Returns:
        float: The Euclidean norm of the wrapped angle vector.
        
    Notes:
        Kinova contains some continuous joints, that's why we need to wrap the angle for some of the joints.
    """
    wrapped_angle = diff
    wrapped_angle[0] = (diff[0] + math.pi) % (2 * math.pi) - math.pi
    wrapped_angle[2] = (diff[2] + math.pi) % (2 * math.pi) - math.pi
    wrapped_angle[4] = (diff[4] + math.pi) % (2 * math.pi) - math.pi
    wrapped_angle[6] = (diff[6] + math.pi) % (2 * math.pi) - math.pi
    return np.linalg.norm(wrapped_angle)

def integrate(model, ts_sim, x0, desired_trajectory, controller, method='RK45'):
    """
    Integrates the dynamics of a robotic system over a given time period.
    Parameters:
    -----------
    model : pinocchio.Model
        The Pinocchio model of the robot, containing its kinematic and dynamic properties.
    ts_sim : array-like
        A sequence of time points at which to solve the dynamics.
    x0 : array-like
        The initial state of the system, consisting of joint positions and velocities.
    desired_trajectory : callable
        A function that takes time `t` as input and returns the desired joint positions,
        velocities, and accelerations as (qd, qd_d, qd_dd).
    controller : callable
        A function that computes the control input (torques) given the current state
        (q, v) and desired trajectory (qd, qd_d, qd_dd).
    method : str, optional
        The numerical integration method to use (default is 'RK45') in scipy.integrate.solve_ivp.
    Returns:
    --------
    qs : ndarray
        The joint positions over time, with shape (len(ts_sim), nq).
    vs : ndarray
        The joint velocities over time, with shape (len(ts_sim), nv).
    taus : ndarray
        The control inputs (torques) over time, with shape (len(ts_sim), nv).
    Notes:
    ------
    - This function uses `scipy.integrate.solve_ivp` to solve the system's dynamics.
    - The dynamics are computed using the Articulated Body Algorithm (ABA) from Pinocchio.
    """
    
    nq = model.nq
    nv = model.nv
    data = model.createData()
    
    # robot dynamics
    def dynamics(t, x):
        """
        Computes the state derivatives for a robotic system using the articulated body algorithm (ABA).
        Parameters:
            t (float): The current time.
            x (numpy.ndarray): The state vector of the system, where the first `nq` elements represent 
                               the joint positions (q) and the remaining elements represent the joint 
                               velocities (v).
        Returns:
            numpy.ndarray: The derivative of the state vector, which includes the joint velocities (v) 
                           and joint accelerations (a).
        """
        
        q = x[:nq]
        v = x[nq:]
        
        qd, qd_d, qd_dd = desired_trajectory(t)
        
        tau = controller(q, v, qd, qd_d, qd_dd)
        
        a = pin.aba(
            model, data, q, v, tau
        )
        
        return np.concatenate([v, a])
    
    # solve the ODE
    print(f"Start integrating using {method} method")
    sol = solve_ivp(
        dynamics, 
        [ts_sim[0], ts_sim[-1]],
        x0, 
        method=method,
        t_eval=ts_sim
    )
    
    # recover the position and velocity
    qs = sol.y[:nq,:].T
    vs = sol.y[nq:,:].T
    
    # recover control input
    taus = np.zeros((len(ts_sim), nv))
    for i, t in enumerate(ts_sim):
        q = qs[i]
        v = vs[i]
        qd, qd_d, qd_dd = desired_trajectory(t)
        taus[i] = controller(q, v, qd, qd_d, qd_dd)
    
    return qs, vs, taus

def set_position(robot, q):
    """
    Sets the joint positions of Kinova in a PyBullet.

    This function iterates through all the joints of the given robot and sets their positions
    based on the provided list of target values `q`. Fixed joints are reset to a target value
    of 0, while other joints are set to the corresponding values in `q`. After updating the
    joint states, the physics simulation is stepped.

    Args:
        robot (int): The unique ID of the robot in the physics simulation.
        q (list of float): A list of target joint positions for the robot. The length of this
            list should match the number of non-fixed joints in the robot.

    Raises:
        IndexError: If the length of `q` is less than the number of non-fixed joints in the robot.
    """
    id = 0
    for i in range(p.getNumJoints(robot)):
        joint_info = p.getJointInfo(robot, i)
        joint_type = joint_info[2]
        if joint_type == p.JOINT_FIXED:
            p.resetJointState(robot, i, targetValue=0)
        else:
            p.resetJointState(robot, i, targetValue=q[id])
            id += 1
    p.stepSimulation()