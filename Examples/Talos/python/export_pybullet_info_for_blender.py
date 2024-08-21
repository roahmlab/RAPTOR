import pickle
import pybullet as p
import numpy as np
import time
from scipy.spatial.transform import Rotation as R

filepath = "../../../Robots/talos/"

# Define the list of meshes with their properties
meshes = [
    {'name': 'base_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/torso/base_link.STL', 'parent': 'Rz', 'material': 'LightGrey'},
    {'name': 'leg_left_1_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/v2/hip_z_lo_res.stl', 'parent': 'leg_left_1_joint', 'material': 'DarkGrey'},
    {'name': 'leg_left_2_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/v2/hip_x_lo_res.stl', 'parent': 'leg_left_2_joint', 'material': 'DarkGrey'},
    {'name': 'leg_left_3_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/v2/hip_y_lo_res.stl', 'parent': 'leg_left_3_joint', 'material': 'DarkGrey'},
    {'name': 'leg_left_4_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/v2/knee_lo_res.stl', 'parent': 'leg_left_4_joint', 'material': 'DarkGrey'},
    {'name': 'leg_left_5_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/v2/ankle_Y_lo_res.stl', 'parent': 'leg_left_5_joint', 'material': 'Grey'},
    {'name': 'leg_left_6_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/v2/ankle_X_lo_res.stl', 'parent': 'leg_left_6_joint', 'material': 'Grey'},
    {'name': 'leg_right_1_link', 'scale': (1.0, -1.0, 1.0), 'file': 'meshes/v2/hip_z_lo_res.stl', 'parent': 'leg_right_1_joint', 'material': 'DarkGrey'},
    {'name': 'leg_right_2_link', 'scale': (1.0, -1.0, 1.0), 'file': 'meshes/v2/hip_x_lo_res.stl', 'parent': 'leg_right_2_joint', 'material': 'DarkGrey'},
    {'name': 'leg_right_3_link', 'scale': (1.0, -1.0, 1.0), 'file': 'meshes/v2/hip_y_lo_res.stl', 'parent': 'leg_right_3_joint', 'material': 'DarkGrey'},
    {'name': 'leg_right_4_link', 'scale': (1.0, -1.0, 1.0), 'file': 'meshes/v2/knee_lo_res.stl', 'parent': 'leg_right_4_joint', 'material': 'DarkGrey'},
    {'name': 'leg_right_5_link', 'scale': (1.0, -1.0, 1.0), 'file': 'meshes/v2/ankle_Y_lo_res.stl', 'parent': 'leg_right_5_joint', 'material': 'Grey'},
    {'name': 'leg_right_6_link', 'scale': (1.0, -1.0, 1.0), 'file': 'meshes/v2/ankle_X_lo_res.stl', 'parent': 'leg_right_6_joint', 'material': 'Grey'},
    {'name': 'torso_1_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/torso/torso_1.STL', 'parent': 'Rz', 'material': 'LightGrey'},
    {'name': 'torso_2_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/torso/torso_2.STL', 'parent': 'Rz', 'material': 'LightGrey'},
    {'name': 'arm_left_1_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/arm/arm_1_collision.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'arm_left_2_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/arm/arm_2_collision.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'arm_left_3_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/arm/arm_3_collision.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'arm_left_4_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/arm/arm_4_collision.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'arm_left_5_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/arm/arm_5_collision.STL', 'parent': 'Rz', 'material': 'LightGrey'},
    {'name': 'arm_left_6_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/arm/arm_6_collision.STL', 'parent': 'Rz', 'material': 'LightGrey'},
    {'name': 'arm_left_7_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/arm/arm_7.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    # {'name': 'wrist_left_ft_link', 'scale': (1.0, 1.0, 1.0), 'file': 'CYLINDER', 'parent': 'Rz', 'material': 'LightGrey'},
    # {'name': 'wrist_left_ft_tool_link', 'scale': (1.0, 1.0, 1.0), 'file': 'CYLINDER', 'parent': 'Rz', 'material': 'DarkerGrey'},
    {'name': 'gripper_left_base_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/base_link.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'gripper_left_inner_double_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/inner_double.STL', 'parent': 'Rz', 'material': 'Orange'},
    {'name': 'gripper_left_fingertip_1_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/fingertip.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'gripper_left_fingertip_2_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/fingertip.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'gripper_left_inner_single_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/inner_single.STL', 'parent': 'Rz', 'material': 'Orange'},
    {'name': 'gripper_left_fingertip_3_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/fingertip.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'gripper_left_motor_double_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/gripper_motor_double.STL', 'parent': 'Rz', 'material': 'Orange'},
    {'name': 'gripper_left_motor_single_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/gripper_motor_single.STL', 'parent': 'Rz', 'material': 'Orange'},
    {'name': 'arm_right_1_link', 'scale': (1.0, -1.0, 1.0), 'file': 'meshes/arm/arm_1_collision.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'arm_right_2_link', 'scale': (1.0, -1.0, 1.0), 'file': 'meshes/arm/arm_2_collision.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'arm_right_3_link', 'scale': (1.0, -1.0, 1.0), 'file': 'meshes/arm/arm_3_collision.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'arm_right_4_link', 'scale': (1.0, -1.0, 1.0), 'file': 'meshes/arm/arm_4_collision.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'arm_right_5_link', 'scale': (1.0, -1.0, 1.0), 'file': 'meshes/arm/arm_5_collision.STL', 'parent': 'Rz', 'material': 'LightGrey'},
    {'name': 'arm_right_6_link', 'scale': (1.0, -1.0, 1.0), 'file': 'meshes/arm/arm_6_collision.STL', 'parent': 'Rz', 'material': 'LightGrey'},
    {'name': 'arm_right_7_link', 'scale': (1.0, -1.0, 1.0), 'file': 'meshes/arm/arm_7_collision.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    # {'name': 'wrist_right_ft_link', 'scale': (1.0, 1.0, 1.0), 'file': 'CYLINDER', 'parent': 'Rz', 'material': 'LightGrey'},
    # {'name': 'wrist_right_ft_tool_link', 'scale': (1.0, 1.0, 1.0), 'file': 'CYLINDER', 'parent': 'Rz', 'material': 'DarkerGrey'},
    {'name': 'gripper_right_base_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/base_link.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'gripper_right_inner_double_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/inner_double.STL', 'parent': 'Rz', 'material': 'Orange'},
    {'name': 'gripper_right_fingertip_1_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/fingertip.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'gripper_right_fingertip_2_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/fingertip.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'gripper_right_inner_single_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/inner_single.STL', 'parent': 'Rz', 'material': 'Orange'},
    {'name': 'gripper_right_fingertip_3_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/fingertip.STL', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'gripper_right_motor_double_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/gripper_motor_double.STL', 'parent': 'Rz', 'material': 'Orange'},
    {'name': 'gripper_right_motor_single_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/gripper/gripper_motor_single.STL', 'parent': 'Rz', 'material': 'Orange'},
    {'name': 'head_1_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/head/head_1.stl', 'parent': 'Rz', 'material': 'DarkGrey'},
    {'name': 'head_2_link', 'scale': (1.0, 1.0, 1.0), 'file': 'meshes/head/head_2.stl', 'parent': 'Rz', 'material': 'LightGrey'},
    # {'name': 'imu_link', 'scale': (1.0, 1.0, 1.0), 'file': 'BOX', 'parent': 'Rz', 'material': 'LightGrey'}
]

# Import robot urdf
# p.connect(p.GUI)
p.connect(p.DIRECT)
urdf_filename = filepath + "talos_reduced_armfixed_floatingbase.urdf"
robot = p.loadURDF(urdf_filename, useFixedBase=True, basePosition=[0,0,0], baseOrientation=[0,0,0,1])
num_joints = p.getNumJoints(robot)
obj_name_to_link_id = {}
obj_name_to_extra_transform = {}

# Import and name each mesh
for mesh in meshes:
    for joint_index in range(num_joints):
        joint_info = p.getJointInfo(robot, joint_index)
        link_name = joint_info[12].decode('utf-8')
        if mesh['name'] == link_name:
            obj_name_to_link_id[mesh['name']] = joint_index
            break
        
    if 'pos' in mesh:
        translation = mesh['pos']
        rotation = R.from_euler('xyz', np.array(mesh['rpy']) * np.pi / 180.0)
        transformation_matrix = np.eye(4)
        transformation_matrix[:3, :3] = rotation.as_matrix()
        transformation_matrix[:3, 3] = translation
        obj_name_to_extra_transform[mesh['name']] = transformation_matrix
        
        
# Read trajectories info
traj_filename = '../data/full-trajectories_forward_0.0.txt'
data = np.loadtxt(traj_filename)
data = data.T

print(data.shape)

# Update joint angles in pybullet and get mesh tranformations
keyframes_obj_name_to_transform = []
for tid in range(0, data.shape[1]):
    pos = data[:18, tid]

    id = 0
    for i in range(num_joints):
        joint_info = p.getJointInfo(robot, i)
        joint_name = joint_info[1].decode('utf-8')
        joint_type = joint_info[2]
        if joint_type == p.JOINT_FIXED:
            if joint_name == 'arm_left_2_joint':
                p.resetJointState(robot, i, targetValue=0.6)
            elif joint_name == 'arm_right_2_joint':
                p.resetJointState(robot, i, targetValue=0.6)
            else:
                p.resetJointState(robot, i, targetValue=0)
        else:
            p.resetJointState(robot, i, targetValue=pos[id])
            id += 1
    p.stepSimulation()
    # time.sleep(5e-2)
            
    # Setup mesh configurations
    obj_name_to_transform = {}
    for obj_name in obj_name_to_link_id:
        if obj_name not in obj_name_to_extra_transform:
            link_id = obj_name_to_link_id[obj_name]
            link_state1 = p.getLinkState(robot, link_id)
            
            obj_location = link_state1[4]
            obj_rotation_quaternion = link_state1[5]
        elif obj_name in obj_name_to_extra_transform:
            link_id = obj_name_to_link_id[obj_name]
            link_state1 = p.getLinkState(robot, link_id)
            location = link_state1[4]
            orientation = link_state1[5]
            
            transformation_matrix = np.eye(4)
            transformation_matrix[:3, :3] = R.from_quat(orientation).as_matrix()
            transformation_matrix[:3, 3] = location
            
            extra_transform = obj_name_to_extra_transform[obj_name]
            
            transformation_matrix = transformation_matrix @ extra_transform
            
            location_new = transformation_matrix[:3, 3]
            orientation_new = R.from_matrix(transformation_matrix[:3, :3]).as_quat()
            
            obj_location = location_new
            obj_rotation_quaternion = orientation_new
            
        obj_name_to_transform[obj_name] = (obj_location, obj_rotation_quaternion)
    
    obj_name_to_transform['Camera'] = ([pos[0] - 2.2, pos[1] - 2.2, 1.8], [np.pi * 0.40, 0, -np.pi*0.25])
    obj_name_to_transform['Area_Light_Up'] = ([pos[0], pos[1], 3], [0, 0, 0])
    # obj_name_to_transform['Area_Light_Forward'] = ([pos[0], pos[1] - 3, 1], [np.pi/2, 0, 0])
        
    keyframes_obj_name_to_transform.append(obj_name_to_transform)
    
# Save the results to a file
with open('../data/keyframes_obj_name_to_transform_forward_0.0.pkl', 'wb') as file:
    pickle.dump(keyframes_obj_name_to_transform, file)
    

# import pinocchio as pin
# from pinocchio.robot_wrapper import RobotWrapper

# def extract_visual_mesh_paths(urdf_file):
#     # Load the robot model from the URDF file
#     robot = RobotWrapper.BuildFromURDF(urdf_file)

#     # Get the visual models
#     visual_model = robot.visual_model

#     # Extract and print mesh file paths
#     mesh_paths = []
#     for geom in visual_model.geometryObjects:
#         scale_str = f"({geom.meshScale[0]}, {geom.meshScale[1]}, {geom.meshScale[2]})"
        
#         if np.abs(geom.meshColor[0] - 0.7) < 1e-6:
#             mesh_color_str = 'Grey'
#         elif np.abs(geom.meshColor[0] - 0.9) < 1e-6:
#             mesh_color_str = 'LightGrey'
#         elif np.abs(geom.meshColor[0] - 0.2) < 1e-6:
#             mesh_color_str = 'DarkGrey'
#         elif np.abs(geom.meshColor[1] - 0.5) < 1e-6:
#             mesh_color_str = 'Orange'
#         else:
#             mesh_color_str = 'DarkerGrey'
        
#         print(f"{{'name': '{geom.name}', 'scale': {scale_str}, 'file': '{geom.meshPath}', 'parent': '{robot.model.names[geom.parentJoint]}', 'material': '{mesh_color_str}'}},")

#     return mesh_paths

# # Example usage:
# mesh_paths = extract_visual_mesh_paths(urdf_filename)