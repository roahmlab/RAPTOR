import sys
sys.path.append("/workspaces/RAPTOR/build/lib")
import end_effector_sysid_nanobind

def main():
    trajectory_filename = "/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_moment_id_data/2024_11_20_gripper_id_20_together.txt"
    urdf_filename  = "/workspaces/RAPTOR/Robots/kinova-gen3/kinova_grasp_fixed.urdf"
    firction_filename = "/workspaces/RAPTOR/Examples/Kinova/SystemIdentification/IterativeSystemIdentification/end_moment_id_data/friction_params.csv"
    print("Testing EndEffectorPybindWrapper")

    wrapper = end_effector_sysid_nanobind.EndEffectorPybindWrapper(
        urdf_filename, trajectory_filename, firction_filename, 10, True
    )
    print("Wrapper initialized successfully.")

    wrapper.set_ipopt_parameters(
        1e-8,            # tol
        100.0,           # max_wall_time
        5,               # print_level
        3000,            # max_iter
        "monotone",      # mu_strategy
        "ma57",          # linear_solver
        False            # gradient_check
    )

    print("IPOPT parameters set successfully.")

    result = wrapper.optimize()
    print("Optimization result:", result)

    distance = wrapper.analyze_solution()
    print("Final distance:", distance)

if __name__ == "__main__":
    main()
