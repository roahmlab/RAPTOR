import pinocchio as pin
import numpy as np
import cvxpy as cp
import sys
import time

# check argue
if len(sys.argv) < 2:
    print("Error: No argument provided. Please specify the file number.")
    sys.exit(1)  

file_number = sys.argv[1]

# Load the robot model
urdf_filename = "/workspaces/RAPTOR/Robots/kinova-gen3/kinova_grasp_fixed.urdf"
model = pin.buildModelFromUrdf(urdf_filename)
data = model.createData()

model.armature.fill(0)
model.damping.fill(0)
model.friction.fill(0)


include_offset_input = False
num_links= model.nv

# load the data
q = np.loadtxt(f"/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_effector_params_data/q_downsampled_{file_number}.csv")
v=  np.loadtxt(f"/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_effector_params_data/q_d_downsampled_{file_number}.csv")
a = np.loadtxt(f"/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_effector_params_data/q_dd_downsampled_{file_number}.csv")
tau_data = np.loadtxt(f"/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_effector_params_data/tau_downsampled_{file_number}.csv")
solution = np.loadtxt("/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/friction_parameters_solution_6.csv")

friction = solution[0 : num_links]
damping = solution [num_links : 2 * num_links]
armature  = solution[2 * num_links : 3 * num_links]

# Compute total friction 
total_friction_force =(
    np.multiply(friction, np.sign(v)) +
    np.multiply(damping, v) +
    np.multiply(armature, a) 
).flatten()

# Compute observation matrix
n_samples = q.shape[0]
Y_data = []
for i in range(n_samples):
    q_i, v_i, a_i = q[i], v[i], a[i]
    Y = pin.computeJointTorqueRegressor(model, data, q_i, v_i, a_i)
    Y_data.append(Y)

# end effector parameters
Y_data = np.vstack(Y_data) 
Y_data_end =Y_data[:,60:70]

tau_data =np.array(tau_data).flatten()

# read urdf inertia parameters 
x = []
for i in range(num_links):
    picc_id = i+1
    x.extend(model.inertias[picc_id].toDynamicParameters())  

z = cp.Variable(Y_data_end.shape[1])  
x_fixed = x[:60]
x_update = cp.hstack([x_fixed, z])


# CVXPY
# Objective function 
objective = cp.Minimize(cp.sum_squares(Y_data @ x_update + total_friction_force - tau_data))

# LMI constraints
constraints = []

# Extract mass, center of mass, and inertia
i = 0
mass = z[10 * i]
com = z[10 * i + 1: 10 * i + 4]
inertia_elements = z[10 * i + 4: 10 * i + 10]

Ixx = inertia_elements[0]
Ixy = inertia_elements[1]
Iyy = inertia_elements[2]
Ixz = inertia_elements[3]
Iyz = inertia_elements[4]
Izz = inertia_elements[5]

inertia = cp.bmat([
    [Ixx, Ixy, Ixz],
    [Ixy, Iyy, Iyz],
    [Ixz, Iyz, Izz]
])


# LMI matrix
LMI = cp.bmat([
    [0.5 * cp.trace(inertia) * np.eye(3) - inertia, cp.reshape(com, (3, 1))],
    [cp.reshape(com, (1, 3)), cp.reshape(mass, (1, 1))]
])

# psd constraints
constraints.append(LMI >> 0)

# solve
start_time = time.time()
prob = cp.Problem(objective, constraints)

prob.solve()
end_time = time.time()
elapsed_time = end_time - start_time
print(f"computation time: {elapsed_time:.2f} s")

# Output results
print("Optimal parameters z:", z.value)
x[60:70] = [round(param, 6) for param in x[60:70]]
print("URDF end parameters:", x[60:70])

# test the result
x_cut = x[0:60]
x_cut.extend(z.value)
x_update = np.array(x_cut)
print("\n",x_update.shape)
print("Y_data @ x_update  + total_friction_force - tau_data:", Y_data @ x_update + total_friction_force- tau_data)

# write to file (test and training same dataset)
outputfile = f"/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_effector_params_data/end_effector_estimate_tau_{file_number}.csv"
np.savetxt(outputfile, (Y_data @ x_update + total_friction_force).reshape(n_samples,7), delimiter=" ")

# write to file (test and training different dataset)
# outputfile1 = f"/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_effector_params_data/end_effector_estimate_tau_{file_number}.csv"
# x_cut = x[0:60]
# x_cut.extend(z.value)
# final = 
# np.savetxt(outputfile1, z.value, delimiter=" ")


# q_test = np.loadtxt(f"/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_effector_params_data/q_downsampled_25.csv")
# v_test=  np.loadtxt(f"/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_effector_params_data/q_d_downsampled_25.csv")
# a_test = np.loadtxt(f"/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_effector_params_data/q_dd_downsampled_25.csv")
# tau_data_test = np.loadtxt(f"/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_effector_params_data/tau_downsampled_25.csv")


# n_samples = q_test.shape[0]
# Y_data_test = []
# for i in range(n_samples):
#     q_i, v_i, a_i = q_test[i], v_test[i], a_test[i]
#     Y = pin.computeJointTorqueRegressor(model, data, q_i, v_i, a_i)
#     Y_data_test.append(Y)
# Y_data_test = np.vstack(Y_data_test) 

# total_friction_force_test =(
#     np.multiply(friction, np.sign(v_test)) +
#     np.multiply(damping, v_test) +
#     np.multiply(armature, a_test) 
# ).flatten()

# outputfile = f"/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_effector_params_data/end_effector_estimate_tau_{file_number}.csv"
# np.savetxt(outputfile, (Y_data_test @ x_update + total_friction_force_test).reshape(n_samples,7), delimiter=" ")

