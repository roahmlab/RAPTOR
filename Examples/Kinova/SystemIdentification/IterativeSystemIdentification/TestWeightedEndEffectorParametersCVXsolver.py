import pinocchio as pin
import numpy as np
import cvxpy as cp
import sys
import time

def optimize_parameters(Y_data_end, weight_tau_end, weigth_extend):
    """
    Performs optimization using CVXPY to estimate the end-effector parameters.

    Parameters:
    - Y_data_end (np.ndarray): Observation matrix for the end-effector parameters.
    - Y_data_fixed (np.ndarray): Observation matrix for the fixed parameters.
    - total_friction_force (np.ndarray): Computed total friction force.
    - tau_data (np.ndarray): Torque data.

    Returns:
    - np.ndarray: Optimized parameters z.
    - float: Time taken to solve the optimization problem.
    """
    # Define CVXPY variable for the end-effector parameters
    z = cp.Variable(Y_data_end.shape[1])

    # Construct the updated parameter vector by concatenating fixed parameters and optimization variables
    objective = cp.Minimize(cp.sum_squares(np.multiply(weigth_extend,Y_data_end) @ z - weight_tau_end))

    # Extract mass, center of mass (com), and inertia elements from z
    mass = z[0]
    com = z[1:4]
    inertia_elements = z[4:10]

    Ixx, Ixy, Iyy, Ixz, Iyz, Izz = inertia_elements

    # Construct the inertia matrix
    inertia = cp.bmat([
        [Ixx, Ixy, Ixz],
        [Ixy, Iyy, Iyz],
        [Ixz, Iyz, Izz]
    ])

    # Construct the Linear Matrix Inequality (LMI)
    LMI = cp.bmat([
        [0.5 * cp.trace(inertia) * np.eye(3) - inertia, cp.reshape(com, (3, 1))],
        [cp.reshape(com, (1, 3)), cp.reshape(mass, (1, 1))]
    ])

    # Define the constraints
    constraints = [LMI >> 0]

    # Define and solve the optimization problem
    prob = cp.Problem(objective, constraints)

    start_time = time.time()
    prob.solve()
    end_time = time.time()

    elapsed_time = end_time - start_time
    print(f"Computation time: {elapsed_time:.2f} s")

    return z.value, elapsed_time

def run_weighted_optimization(n_samples, num_links, x_fixed , weight_tol, o_sqrt_tol, tau_data, Y_data, total_friction_force):
    weighted_old = np.ones(n_samples *num_links)*10000
    weighted_new = np.ones(n_samples*num_links)
    weigth_extend = np.tile(weighted_new[:, np.newaxis], (1, 10))
   

    error_threshold = 5

    x_sol = np.zeros(10)
    Rh = np.zeros(n_samples * num_links)

    while np.linalg.norm(weighted_new - weighted_old) > weight_tol:
        print("weight differece", np.linalg.norm(weighted_new - weighted_old))
        Omega_old = np.eye(num_links)*10000
        Omega_new = np.eye(num_links)
        Omega_sqrt = np.linalg.cholesky(Omega_new)

        while(np.linalg.norm(Omega_new - Omega_old) > o_sqrt_tol):
            # print("differece", np.linalg.norm(Omega_new - Omega_old))
            Ts = np.zeros(n_samples* num_links)
            friction_s = np.zeros(n_samples * num_links)
            Ys = np.zeros((n_samples* num_links, Y_data.shape[1]))
            for i in range(n_samples):
                Ts[i * num_links:(i + 1) * num_links] = (
                                    np.linalg.solve(Omega_sqrt, tau_data[i * num_links:(i + 1) * num_links]))

                friction_s[i * num_links:(i + 1) * num_links] = (
                                    np.linalg.solve(Omega_sqrt,  total_friction_force[i * num_links:(i + 1) * num_links]))

                Ys[i * num_links:(i + 1) * num_links, :] =( 
                                    np.linalg.solve(Omega_sqrt, Y_data[i * num_links:(i + 1) * num_links, :]))


            # Torque and regressor after weight
            tau_end = Ts- (Ys[:,:60] @ x_fixed + friction_s) 
            weight_tau_end = np.multiply(weighted_new,tau_end)
         
            z_opt, _= optimize_parameters(Ys[:,60:70], weight_tau_end, weigth_extend)
            # print("z_opt", z_opt)

            x_sol = z_opt
            Rh = weight_tau_end - np.multiply(weigth_extend,Ys[:,60:70]) @ x_sol
            # print("Rh", Rh)
            Eh = Rh.reshape((num_links, n_samples))
            # print("Eh", Eh)


            # Update covariance matrix
            Omega_old = Omega_new
            # Omega_new = (Omega_sqrt @ Eh @ Eh.T @ Omega_sqrt) / (n_samples- 10)
            Omega_new = (Eh @ Eh.T) / (n_samples - 10)
            # print("Omega_new", Omega_new)
            epsilon = 1e-6  # Regularization parameter
            Omega_new += epsilon * np.eye(num_links)

            # Update sqrt covariance matrix
            Omega_sqrt = np.linalg.cholesky(Omega_new)

            if not np.all(np.isfinite(Omega_sqrt)):
                break

        # Update weight
        weighted_old = weighted_new
        Psi = np.where(np.abs(Rh) <= error_threshold, 1, 0)
        weighted_new = np.minimum(weighted_old, Psi)
        print("weighted_new", weighted_new)

        # Update weight_extend
        weight_extend = np.tile(weighted_new[:, np.newaxis], (1, 10))

    return x_sol



def main():
    """
    Main function to execute the system identification process.
    """
    # Check for the required command-line argument
    if len(sys.argv) < 2:
        print("Error: No argument provided. Please specify the file number.")
        sys.exit(1)  
    
    file_number = sys.argv[1]

    # Load the robot model from URDF
    urdf_filename = "/workspaces/RAPTOR/Robots/kinova-gen3/kinova_grasp_fixed.urdf"
    model = pin.buildModelFromUrdf(urdf_filename)
    data = model.createData()

    # Initialize model parameters
    model.armature.fill(0)
    model.damping.fill(0)
    model.friction.fill(0)

    num_links = model.nv  # Number of velocity degrees of freedom

    # Load data from CSV files
    base_path = "/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/"
    q = np.loadtxt(f"{base_path}end_effector_params_data/q_downsampled_{file_number}.csv")
    v = np.loadtxt(f"{base_path}end_effector_params_data/q_d_downsampled_{file_number}.csv")
    a = np.loadtxt(f"{base_path}end_effector_params_data/q_dd_downsampled_{file_number}.csv")
    tau_data = np.loadtxt(f"{base_path}end_effector_params_data/tau_downsampled_{file_number}.csv")
    solution = np.loadtxt(f"{base_path}full_params_data/friction_parameters_solution_6.csv")

    tau_data = tau_data.flatten()
    # Extract friction, damping, and armature parameters from the solution
    friction = solution[0 : num_links]
    damping = solution[num_links : 2 * num_links]
    armature = solution[2 * num_links : 3 * num_links]

    # Compute total friction force
    total_friction_force = (
        np.multiply(friction, np.sign(v)) +
        np.multiply(damping, v) +
        np.multiply(armature, a)
    ).flatten()

    # Compute the observation matrix Y_data
    n_samples = q.shape[0]
    Y_data = []
    for i in range(n_samples):
        q_i, v_i, a_i = q[i], v[i], a[i]
        Y = pin.computeJointTorqueRegressor(model, data, q_i, v_i, a_i)
        Y_data.append(Y)

    # Load URDF inertia parameters
    x = []
    for i in range(num_links):
        picc_id = i + 1
        x.extend(model.inertias[picc_id].toDynamicParameters())  

    # Separate the observation matrix for end-effector parameters and fixed parameters
    Y_data = np.vstack(Y_data) 
    x_fixed = np.array(x[:60])


    # Y_data_end = Y_data[:,60:70]
    # Y_data_fixed = Y_data[:, :60]

    # tau_fixed = Y_data_fixed @ x_fixed
    weight_tol, o_sqrt_tol = 4e-0, 4e-0


    # Perform optimization to estimate end-effector parameters
    # z_opt, computation_time = optimize_parameters(Y_data_end, tau_fixed, total_friction_force, tau_data)
    z_opt = run_weighted_optimization(n_samples, num_links, x_fixed , weight_tol, o_sqrt_tol, tau_data, Y_data, total_friction_force)
    # Output the optimized parameters
    print("Optimal parameters z:", z_opt)
    urdf_end_params = np.round(x[60:70], 6)
    print("URDF end parameters:", urdf_end_params)

    # Update the full parameter vector with the optimized z
    # x_update = np.concatenate([x_fixed, z_opt])
    # print("\nParameter vector shape:", x_update.shape)
    # residual = Y_data @ x_update + total_friction_force - tau_data
    # print("Residual (Y_data @ x_update + total_friction_force - tau_data):", residual)

    # Save the estimated torque to a CSV file
    # outputfile = f"{base_path}end_effector_estimate_tau_{file_number}.csv"
    # np.savetxt(outputfile, (Y_data @ x_update + total_friction_force).reshape(n_samples,7), delimiter=" ")
    # print(f"Estimated torque saved to {outputfile}")

    # Optionally, save the optimized z parameters separately
    # Uncomment the following lines if needed
    # outputfile_z = f"{base_path}end_effector_parameters_z_{file_number}.csv"
    # np.savetxt(outputfile_z, z_opt, delimiter=" ")
    # print(f"Optimized z parameters saved to {outputfile_z}")

if __name__ == "__main__":
    main()


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

