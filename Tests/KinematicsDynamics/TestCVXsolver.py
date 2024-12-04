import pinocchio as pin
import numpy as np
import cvxpy as cp

# Load the robot model
urdf_filename = "/workspaces/RAPTOR/Robots/kinova-gen3/kinova.urdf"
model = pin.buildModelFromUrdf(urdf_filename)
data = model.createData()

model.armature.fill(0)
model.damping.fill(0)
model.friction.fill(0)

# Randomly generate 500 
n_samples = 500
num_links= model.nv
print(num_links)

q = 2*np.pi*  (np.random.rand(n_samples, num_links)-0.5)
v = 2*np.pi* (np.random.rand(n_samples, num_links)-0.5)
a = 2*np.pi* (np.random.rand(n_samples, num_links)-0.5)
# # Store observation matrices and torques
Y_data = []
tau_data = []

for i in range(n_samples):
    q_i, v_i, a_i = q[i], v[i], a[i]
    # Compute observation matrix
    Y = pin.computeJointTorqueRegressor(model, data, q_i, v_i, a_i)
    Y_data.append(Y)
    
    # Compute torque
    tau = pin.rnea(model, data, q_i, v_i, a_i)
    tau_data.append(tau)

# full parameters
Y_data = np.vstack(Y_data) 
tau_data = np.array(tau_data).flatten()


x = []
for i in range(num_links):
    picc_id = i+1
    x.extend(model.inertias[picc_id].toDynamicParameters())          


z = cp.Variable(Y_data.shape[1])  
print(Y_data.shape[1])


objective = cp.Minimize(cp.sum_squares(Y_data @ z - tau_data))

# LMI constraints
constraints = []

for i in range(num_links):
    # Extract mass, center of mass, and inertia
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

    # psd
    constraints.append(LMI[0, 0] >= 0) #scalar
    constraints.append(LMI[:1, :1] >> 0)
    constraints.append(LMI[:2, :2] >> 0)
    constraints.append(LMI>> 0)

# solve
prob = cp.Problem(objective, constraints)
prob.solve()


# Output results
print("Optimal parameters z:", z.value)
print("Urdf parameters x: ", x)
print("Y_data @ z - tau_data:", Y_data @ z.value- tau_data)
