import pickle
import pybullet as p
import numpy as np
import time
from scipy.spatial.transform import Rotation as R

filepath = "../../../Robots/digit-v3/"

# Define the list of meshes with their properties
meshes = [
    {'name': 'torso/base-0', 'scale': (1, 1, 1), 'file': 'meshes/torso-v3/torso.stl', 'parent': 'Rz', 'pos': [0, 0, 0], 'rpy': [0, 0, 0]},
    {'name': 'torso/base-1', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/hip-roll-housing.stl', 'parent': 'Rz', 'pos': [-0.001, 0.091, 0], 'rpy': [0, -90, 0]},
    {'name': 'torso/base-2', 'scale': (1, 1, 1), 'file': 'meshes/arm-v3/link-0.stl', 'parent': 'Rz', 'pos': [-0.001, 0.12, 0.4], 'rpy': [-9.999933098, -90, 0]},
    {'name': 'torso/base-3', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/hip-roll-housing.stl', 'parent': 'Rz', 'pos': [-0.001, -0.091, 0], 'rpy': [0, -90, 0]},
    {'name': 'torso/base-4', 'scale': (1, -1, 1), 'file': 'meshes/arm-v3/link-0.stl', 'parent': 'Rz', 'pos': [-0.001, -0.12, 0.4], 'rpy': [9.999933098, -90, 0]},
    {'name': 'left-leg/hip-roll-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/hip-yaw-housing.stl', 'parent': 'left_hip_roll'},
    {'name': 'left-leg/hip-yaw-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/hip-pitch-housing.stl', 'parent': 'left_hip_yaw'},
    {'name': 'left-leg/hip-pitch-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/hip-pitch.stl', 'parent': 'left_hip_pitch'},
    {'name': 'left-leg/achilles-rod-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/achilles-rod.stl', 'parent': 'left_achilles_rod'},
    {'name': 'left-leg/knee-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/knee.stl', 'parent': 'left_knee'},
    {'name': 'left-leg/shin-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/shin.stl', 'parent': 'left_shin'},
    {'name': 'left-leg/tarsus-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/tarsus.stl', 'parent': 'left_tarsus'},
    {'name': 'left-leg/heel-spring-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/heel-spring.stl', 'parent': 'left_heel_spring'},
    {'name': 'left-leg/toe-a-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/toe-output.stl', 'parent': 'left_toe_A'},
    {'name': 'left-leg/toe-a-rod-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/toe-a-rod.stl', 'parent': 'left_toe_A_rod'},
    {'name': 'left-leg/toe-b-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/toe-output.stl', 'parent': 'left_toe_B'},
    {'name': 'left-leg/toe-b-rod-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/toe-b-rod.stl', 'parent': 'left_toe_B_rod'},
    {'name': 'left-leg/toe-pitch-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/toe-pitch.stl', 'parent': 'left_toe_pitch'},
    {'name': 'left-leg/toe-roll-0', 'scale': (1, 1, 1), 'file': 'meshes/leg-v3/toe-roll.stl', 'parent': 'left_toe_roll'},
    {'name': 'left-arm/shoulder-roll-0', 'scale': (1, 1, 1), 'file': 'meshes/arm-v3/link-1.stl', 'parent': 'left_shoulder_roll'},
    {'name': 'left-arm/shoulder-pitch-0', 'scale': (1, 1, 1), 'file': 'meshes/arm-v3/link-2.stl', 'parent': 'left_shoulder_pitch'},
    {'name': 'left-arm/shoulder-yaw-0', 'scale': (1, 1, 1), 'file': 'meshes/arm-v3/link-3.stl', 'parent': 'left_shoulder_yaw'},
    {'name': 'left-arm/elbow-0', 'scale': (1, 1, 1), 'file': 'meshes/arm-v3/link-4.stl', 'parent': 'left_elbow'},
    {'name': 'right-leg/hip-roll-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/hip-yaw-housing.stl', 'parent': 'right_hip_roll'},
    {'name': 'right-leg/hip-yaw-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/hip-pitch-housing.stl', 'parent': 'right_hip_yaw'},
    {'name': 'right-leg/hip-pitch-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/hip-pitch.stl', 'parent': 'right_hip_pitch'},
    {'name': 'right-leg/achilles-rod-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/achilles-rod.stl', 'parent': 'right_achilles_rod'},
    {'name': 'right-leg/knee-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/knee.stl', 'parent': 'right_knee'},
    {'name': 'right-leg/shin-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/shin.stl', 'parent': 'right_shin'},
    {'name': 'right-leg/tarsus-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/tarsus.stl', 'parent': 'right_tarsus'},
    {'name': 'right-leg/heel-spring-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/heel-spring.stl', 'parent': 'right_heel_spring'},
    {'name': 'right-leg/toe-a-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/toe-output.stl', 'parent': 'right_toe_A'},
    {'name': 'right-leg/toe-a-rod-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/toe-a-rod.stl', 'parent': 'right_toe_A_rod'},
    {'name': 'right-leg/toe-b-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/toe-output.stl', 'parent': 'right_toe_B'},
    {'name': 'right-leg/toe-b-rod-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/toe-b-rod.stl', 'parent': 'right_toe_B_rod'},
    {'name': 'right-leg/toe-pitch-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/toe-pitch.stl', 'parent': 'right_toe_pitch'},
    {'name': 'right-leg/toe-roll-0', 'scale': (1, -1, 1), 'file': 'meshes/leg-v3/toe-roll.stl', 'parent': 'right_toe_roll'},
    {'name': 'right-arm/shoulder-roll-0', 'scale': (1, -1, 1), 'file': 'meshes/arm-v3/link-1.stl', 'parent': 'right_shoulder_roll'},
    {'name': 'right-arm/shoulder-pitch-0', 'scale': (1, -1, 1), 'file': 'meshes/arm-v3/link-2.stl', 'parent': 'right_shoulder_pitch'},
    {'name': 'right-arm/shoulder-yaw-0', 'scale': (1, -1, 1), 'file': 'meshes/arm-v3/link-3.stl', 'parent': 'right_shoulder_yaw'},
    {'name': 'right-arm/elbow-0', 'scale': (1, -1, 1), 'file': 'meshes/arm-v3/link-4.stl', 'parent': 'right_elbow'}
]

# Import robot urdf
# p.connect(p.GUI)
p.connect(p.DIRECT)
urdf_filename = filepath + "digit-v3-armfixedspecific-floatingbase-springfixed.urdf"
robot = p.loadURDF(urdf_filename, useFixedBase=True, basePosition=[0,0,0], baseOrientation=[0,0,0,1])
num_joints = p.getNumJoints(robot)
obj_name_to_link_id = {}
obj_name_to_extra_transform = {}

# Import and name each mesh
for mesh in meshes:
    for joint_index in range(num_joints):
        joint_info = p.getJointInfo(robot, joint_index)
        link_name = joint_info[12].decode('utf-8')
        if mesh['parent'] == link_name:
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
traj_filename = '../data/full-trajectories_forward_1.0.txt'
data = np.loadtxt(traj_filename)
data = data.T

print(data.shape)

# Update joint angles in pybullet and get mesh tranformations
keyframes_obj_name_to_transform = []
for tid in range(0, data.shape[1]):
    pos = data[:36, tid]

    id = 0
    for i in range(num_joints):
        joint_info = p.getJointInfo(robot, i)
        joint_type = joint_info[2]
        if joint_type == p.JOINT_FIXED:
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
        
    obj_name_to_transform['Camera'] = ([pos[0] - 2.2, pos[1] - 2.2, 1.6], [np.pi * 0.40, 0, -np.pi*0.25])
    obj_name_to_transform['Area_Light_Up'] = ([pos[0], pos[1], 3], [0, 0, 0])
    obj_name_to_transform['Area_Light_Forward'] = ([pos[0], pos[1] - 3, 1], [np.pi/2, 0, 0])
    
    # obj_name_to_transform['Camera'] = ([-2.2, -2.2, 1.6], [np.pi * 0.42, 0, -np.pi*0.25])
    # obj_name_to_transform['Area_Light_Up'] = ([0, 0, 3], [0, 0, 0])
    # obj_name_to_transform['Area_Light_Forward'] = ([0, -3, 1], [np.pi/2, 0, 0])
        
    keyframes_obj_name_to_transform.append(obj_name_to_transform)
    
# Save the results to a file
with open('../data/keyframes_obj_name_to_transform_upstairs.pkl', 'wb') as file:
    pickle.dump(keyframes_obj_name_to_transform, file)
    

