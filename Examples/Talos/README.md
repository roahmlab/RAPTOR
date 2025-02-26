# Gait Optimization For Talos

This folder contains the gait optimization implementation for Talos, that is used for comparison between RAPTOR and [aligator](https://github.com/Simple-Robotics/aligator/).

## Getting Started

Our main function is mainly at `TalosSingleStepFixedPosition.cpp`, where we try to solve for a trajectory that makes Talos walk forward with a user-specified step length while starting from a fixed initial configuration.

Create a `data/` folder first in `Examples/Talos` to store all the results.

Specify the step length in `singlestep_optimization_settings.yaml` and run the following command in `build/`:
```bash
./Talos_example
```
This program runs the optimization and then save the trajectories and the corresponding inverse dynamics control inputs in `data/solution-talos-forward-[step length].txt`, over a much finer time discretization.
Check line 196-210 in `TalosSingleStepFixedPosition.cpp`.

In `Examples/Talos/python/talos_simulation.py`, you will be able to find how we simulate the Talos using `scipy.solve_ivp` as the integrator and pinocchio to evaluate the contact dynamics.
Here we apply a tracking controller in the simulation:
```math
u(t) = u_{openloop}(t) + K_Pe(t) + K_D\dot{e}(t)
```
Here, $u_{openloop}(t)$ represents the open loop control input recovered from the control input sequence saved in `data/solution-talos-forward-[step length].txt` using zero-order hold ([ZOH](https://en.wikipedia.org/wiki/Zero-order_hold)).
We also added a simple PD control here where the desired trajectories are also recovered from `data/solution-talos-forward-[step length].txt` using ZOH.

For our implementation on the comparisons using aligator, please refer to our [fork](https://github.com/roahmlab/aligator-roahmlab/) on aligator.

## Migrating to Other Humanoid Robots

This folder provides a framework for optimizing gaits for different humanoid robots. To migrate the code, you can copy this folder and replace all instances of `Talos` with the name of your robot.

Below are the key modifications required for migration:

---

### **URDF Modifications**

1. **URDF File Paths**  
   Update `urdf_filename` in the main files such as `TalosSingleStep.cpp` and `TalosMultipleStep.cpp`.  
   You may also need to update the `filepath` in the same files.

2. **Floating Base Representation**  
   RAPTOR requires a URDF where the robot is attached to a **floating base** with:
   - Three **prismatic joints** (for translation)
   - Three **revolute joints** (for rotation)

   You need to **manually create** this modified URDF.  
   Refer to [talos_reduced_armfixed_floatingbase.urdf](../../Robots/talos/talos_reduced_armfixed_floatingbase.urdf) as an example.

3. **Torso Rotation Joint Name**  
   In **line 214** of `TalosDynamicsConstraints.cpp`, ensure that the `Rz` joint (which allows torso rotation around the **world z-axis**) matches the corresponding joint in your URDF.

   - For example, in the same location of [`G1DynamicsConstraints.cpp`](../Unitree-G1/src/G1DynamicsConstraints.cpp), this joint is renamed to `Rz_joint` to be consistent with the corresponding [URDF](../../Robots/unitree-g1/g1_22dof_floatingbase.urdf).

---

### **Configuration Files**
Update the paths to the optimization configuration files in the main files:

- `singlestep_optimization_settings.yaml`
- `multiplestep_optimization_settings.yaml`

These paths appear in `TalosSingleStep.cpp` and `TalosMultipleStep.cpp`.

---

### **Robot-Related Constants**
Modify the following constants in `TalosConstants.h` to match your robot:

1. **Joint Limits & Torque Limits**  
   These should be available from your robot's URDF.

2. **Ground Friction Coefficients**  
   - `MU`: Translational friction coefficient. 
   - `GAMMA`: Torsional friction coefficient.

---

### **Feet (Ground Contact Frames)**

1. **Foot Joint Names**  
   Update `LEFT_FOOT_NAME` and `RIGHT_FOOT_NAME` in `TalosConstants.h`.  
   These are typically the last joints in the kinematic chain of each leg.

2. **Foot Frame Transformation**  
   The **foot frame** refers to the **center of the bottom of the foot**, which may not coincide with the foot joint.  
   Update the transformation matrices in:
   - `stance_foot_endT` in `TalosDynamicsConstraints.cpp`
   - `leftfoot_endT` and `rightfoot_endT` in `TalosCustomizedConstraints.cpp`

   **Where to find this transformation matrix?**  
   - **Digit-v3**: Available in its MuJoCo XML file (`digit-v3.xml`).  
   - **Talos**: Defined in the URDF (same as the foot link mesh frame).  
   - **Unitree-G1**: Defined in the URDF (foot frame is at the center of four corner spheres used for collision).  

3. **Foot Size**  
   Update `FOOT_WIDTH` and `FOOT_LENGTH` in `TalosConstants.h`.  
   - `FOOT_WIDTH` = **half** of the size of foot along the x axis of the foot frame.
   - `FOOT_LENGTH` = **half** of the size of foot along the y axis of the foot frame.

   Note that `FOOT_WIDTH` is not necessarily the short side and `FOOT_LENGTH` is not necessarily the long side.
   If the robot exhibits instability (e.g., tilting unnaturally), you may have assigned incorrect foot dimensions.

---

By making these modifications, you should be able to optimize gaits for your humanoid robot using this framework.
One example of this migration can be found in [Unitree-G1](../Unitree-G1/), where we duplicated the contents of this folder and modified the necessary components as outlined above.
