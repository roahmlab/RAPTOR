import numpy as np
import sys
import matplotlib.pyplot as plt
import time
import pybullet as p

from kinova_dynamics import goal_distance, set_position

sys.path.append("/workspaces/RAPTOR/build/lib")
import armour_nanobind as armour
import KinovaHLP_nanobind as HLP
import end_effector_sysid_momentum_nanobind

### initializations
display_info = False

urdf_filename = "/workspaces/RAPTOR/Robots/kinova-gen3/gen3_2f85_fixed_object.urdf" # this urdf only has 7 joints for visualization
config_filename = "/workspaces/RAPTOR/Examples/Kinova/Armour/KinovaWithGripperInfo.yaml" # edit this yaml file if necessary

p.connect(p.GUI)
robot = p.loadURDF(urdf_filename, useFixedBase=True)
num_joints = p.getNumJoints(robot)

rrtplanner = HLP.WaypointPlanningPybindWrapper(urdf_filename)
planner = armour.ArmourPybindWrapper(urdf_filename, config_filename, display_info)

# obstacle information (xyz, rpy, size)
obstacles = np.array([[0.4, 0.0, 0.15, 0.0, 0.0, 0.0, 0.30, 0.30, 0.30], # another obstacle in the middle
                      [0.4, 0.0, 0.45, 0.0, 0.0, 0.0, 0.30, 0.30, 0.30], # another obstacle in the middle
                      [0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.01], # ground
                      [-0.39, -0.2, 0.56, 0.0, 0.0, 0.0, 0.08, 2.0, 1.12], # side wall
                     ])

# start and goal in joint configuration space
start = np.array([-1.16396061,  0.34987045, -3.68094828, -1.75466411,  0.19968747, -1.09177839, -0.19841415])
goal = np.array([-0.66290119, -0.09099118, -3.48809793, -0.79306897, 0.04506539, -1.45191409, 0.6983348])
# goal = np.array([-0.02500949, 0.55020233, -3.05694277, -0.54403656, -0.04979576, -2.04844549, 1.59508374])

# initialize the robot to start configuration
set_position(robot, start)

# create and visualize the obstacles in pybullet
for i, obs in enumerate(obstacles):
    position = obs[0:3]
    orientation = p.getQuaternionFromEuler(obs[3:6])
    size = obs[6:9]
    
    boxCollisionShapeId = p.createCollisionShape(shapeType=p.GEOM_BOX,
                                                 halfExtents=size*0.5)
    boxVisualShapeId = p.createVisualShape(shapeType=p.GEOM_BOX,
                                           halfExtents=size*0.5,
                                           rgbaColor=[1,0,0,0.2])
    box_id = p.createMultiBody(baseMass=0,
                               baseVisualShapeIndex=boxVisualShapeId,
                               baseCollisionShapeIndex=boxCollisionShapeId,
                               basePosition=position,
                               baseOrientation=orientation)

### initialization for system identification 
# initial information of the end effector inertial parameters 
theta_lb = np.array([1.2, -0.4, -0.4, -1.0, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2])
theta_ub = np.array([5.2,  0.4,  0.4, -0.1,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2])
                   # mass, comx, comy, comz, Ixx,  Ixy,  Iyy,  Ixz,  Iyz,  Izz
theta_estimation = 0.5 * (theta_lb + theta_ub)
    
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
    10,                  # forward integration horizon
    'Second',            # format of time data, 'Second' or 'Nanosecond'
    True                 # display information or not
)
sysid_solver.set_ipopt_parameters(
    1e-14,          # tol
    5.0,            # max_wall_time
    5,              # print_level
    500,            # max_iter
    "adaptive",     # mu_strategy
    "ma86",         # linear_solver
    False           # gradient_check
)

### start local planning
duration = 3.0 # sec, duration of the (Bezier) trajectory
dt = 1e-4 # sec, sampling frequency

# set up local planner (dump all the information into it)
planner.set_obstacles(obstacles)
planner.set_ipopt_parameters(
    1e-6,           # tol
    1e-6,           # constr_viol_tol
    10.0,           # obj_scaling_factor
    5.0,            # max_wall_time
    5,              # print_level
    "adaptive",     # mu_strategy
    "ma57",         # linear_solver
    False           # gradient_check
)    

# initial conditions of the robot   
q0 = start
qd0 = np.zeros(7)
qdd0 = np.zeros(7)
 
# choose to minimize the distance between the goal and which part of the trajectory
# for example, 0.0 means minimizing the distance between the goal and the start of the trajectory
# 1.0 means minimizing the distance between the goal and the end of the trajectory
t_plan = 1.0 # choose between [0.0, 1.0]

for iter in range(20):
    # set up the local planner
    q_des = goal # desired joint configuration of the current planning iteration
    k_center = np.zeros(7)
    k_range = np.ones(7)
    k_range[:4] = np.pi/64
    k_range[4:] = np.pi/24
    
    planner.set_trajectory_parameters(
        q0, \
        qd0, \
        qdd0, \
        k_center, \
        k_range, \
        duration, \
        q_des, \
        t_plan
    )
    planner.set_endeffector_inertial_parameters(
        theta_estimation,
        theta_lb,
        theta_ub
    )

    # optimize the trajectory
    k, ifFeasible = planner.optimize()
    print(k)

    if ifFeasible:
        info, \
        sphere_x, sphere_y, sphere_z, sphere_radii, \
        torque_centers, torque_radii = \
            planner.analyze_solution()
        
        sphere_info = np.stack((sphere_x, sphere_y, sphere_z), axis=0)
        previous_info = info
            
        # visualize the trajectory
        for tid in range(info.shape[1]):
            set_position(robot, info[:7, tid])
            time.sleep(0.01)
    else:
        print("Failed to find a feasible trajectory")
        break
    
    # print trajectory data in file
    planner.get_trajectory_data(
        np.arange(0, duration, dt),
        "traj_data.csv"
    )
    
    # solve system identification problem by adding the latest trajectory data
    sysid_solver.add_trajectory_file("traj_data.csv")
    theta_solution, theta_uncertainty = sysid_solver.optimize()
    
    print("end effector inertial parameters:")
    print("solution:\n", theta_solution)
    print("uncertainty:\n", theta_uncertainty)
    
    theta_lb_new = theta_solution - theta_uncertainty
    theta_ub_new = theta_solution + theta_uncertainty
    
    for i in range(10):
        if theta_lb_new[i] > theta_lb[i]:
            theta_lb[i] = theta_lb_new[i]
        if theta_ub_new[i] < theta_ub[i]:
            theta_ub[i] = theta_ub_new[i]
        theta_estimation[i] = 0.5 * (theta_lb[i] + theta_ub[i])
        
    print("lower bound:\n", theta_lb)
    print("upper bound:\n", theta_ub)
    
    # update robot information for the next planning iteration
    q0 = info[:7, -1]
    qd0 = np.zeros(7)
    qdd0 = np.zeros(7)
    
    if goal_distance(q0 - goal) < 0.1:
        print("Goal reached")
        break
        
p.disconnect()