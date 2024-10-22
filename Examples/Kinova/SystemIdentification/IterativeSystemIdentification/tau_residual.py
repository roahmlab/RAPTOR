import numpy as np
import matplotlib.pyplot as plt

# Load data
num_joints = 7
file_t = 'full_params_data/t_downsampled_5.csv'
# file_estimate_tau = 'full_params_data/full_parameters_estimate_tau_5.csv'
file_estimate_tau = 'full_params_data/friction_estimate_tau_5.csv'

file_ground_tau = 'full_params_data/tau_downsampled_5.csv'

# Read data from CSV files
t = np.loadtxt(file_t, delimiter=' ')
estimate_tau = np.loadtxt(file_estimate_tau, delimiter=' ')
tau_ground = np.loadtxt(file_ground_tau, delimiter=' ')
tau_ground = -tau_ground

# Create a figure
fig, axs = plt.subplots(3, 3, figsize=(15, 10))
fig.suptitle('data 5 Estimation vs. Ground Truth', fontsize=16)

# Iterate through each joint to create subplots
for k in range(num_joints):
    row = k // 3
    col = k % 3
    ax = axs[row, col]
    
    ax.plot(t, estimate_tau[:, k], "r", label='Estimated')
    ax.plot(t, tau_ground[:, k], "b", label='Ground Truth')
    
    ax.set_title(f'Joint #{k + 1}')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Torque (Nm)')
    ax.legend(loc='best')

# Hide any unused subplots (if less than 9)
for i in range(num_joints, 9):
    fig.delaxes(axs.flat[i])

# Adjust layout and show plot
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
