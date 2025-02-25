# Persistently Exciting Trajectories for System Identification

## Introduction

Persistently exciting trajectories are crucial in robot system identification because they ensure that the collected data sufficiently excites all the relevant system dynamics, allowing for accurate parameter estimation. 
Without persistent excitation, the system identification process may result in ill-conditioned models that fail to capture the full behavior of the robot [1].
The most common way to generate exciting trajectories is to formulate an optimization problem to minimize the condition number of the regressor matrix [2].

## Highlights in RAPTOR

### Analytical gradient for faster deployment
The gradient of such an optimization problem is typically derived using symbolic expressions in MATLAB or Mathematica, making it difficult to generalize across different robots. 
[3] is the first work to derive the analytical gradient of the condition number of the regressor matrix. 
**RAPTOR** implements this first-order analytical gradient in the `RegressorInverseDynamics` class, allowing the code to seamlessly adapt to different robots as long as a URDF is provided.

### Different trajectory primitives
As introduced in the `Trajectories` class, various trajectory primitives can be used to generate exciting trajectories.
The most common choice is `FixedFrequencyFourierCurves`.
An alternative option could be `PiecewiseBezierCurves` as introduced in [4].

### Collision Avoidance
Very few related works consider collision avoidance when designing exciting trajectories, even though it may be unavoidable in certain settings. 
To the best of our knowledge, [6] is the only work that introduces collision avoidance in trajectory excitation; however, it focuses solely on path planning rather than full trajectory planning. 
**RAPTOR** integrates collision avoidance constraints directly into the optimization problem, enhancing safety across different workspace settings.

## Getting Started

There are two examples: `KinovaRegressorExample.cpp` and `KinovaRegressorExamplePayload.cpp`, each designed for different applications.

In the trajectory optimization problem, the following constraints are enforced, with first-order analytical gradients provided:
 - TrajectoryTerminalConstraints: ensures zero velocity and acceleration at the end for safety
 - Joint limits
 - Velocity limits
 - Acceleration limits
 - Torque limits (based on the URDF from the manufacturer)
 - Collision avoidance with box obstacles

We choose the `FixedFrequencyFourierCurves` class, which also allows users to specify the initial position and velocity as an option. 
The initial velocity is typically set to zero, as the robot generally starts its motion from rest.

### Whole-body Exciting Trajectory Generation
`KinovaRegressorExample.cpp` implements the trajectory optimization problem that minimizes the condition number of the regressor matrix, used to identify the inertial parameters of the whole robot.

Note that the regressor matrix is typically **rank-deficient** for nearly all robots, meaning some inertial parameters cannot be independently identified. 
This is where the concept of **base inertial parameters** comes into play [2].
As a result, we only need to optimize the condition number of the **linearly independent** part of the regressor matrix.
In this case, the matrix is full rank, ensuring that the condition number remains finite. 
The `independent_param_inds` variable defined in `KinovaRegressorExample.cpp` specifies the index of the linearly independent columns of the regressor matrix, which also corresponds to the identifiable inertial parameters of the robot.

To identify the base inertial parameters of the robot, we strongly recommend that users read [7], which provides a comprehensive overview of the methodology.
The authors also provide an open-source implementation at [spatial_v2_extended](https://github.com/ROAM-Lab-ND/spatial_v2_extended/).
Users should refer to this script [RPNA_URDF_Example.m](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/identifiability/RPNA_URDF_Example.m).
The variable `ind` stores the indices of identifiable inertial parameters in each link.
Based on this example script, we would suggest the following changes at the end:
```MATLAB
sympref('FloatingPointOutput',true);
independent_param_inds = [];
for i = 1:size(Perp_Basis_sym,2)
    ind = find(Perp_Basis_sym(:,i)==1,1);
    independent_param_inds = [independent_param_inds, ind];
    
    % Identifable parameter combination
    sym_result = Perp_Basis_sym(:,i)'*sym_params;
    
    % Work to strip out zero coefficients
    [coef, monomials] = coeffs(sym_result);
    coef = CleanMat(coef);
    sym_result = simplify( sum( coef(:).*monomials(:) ) );
    
    % And then group terms that multiply each parameter
    sym_result = jacobian(sym_result, sym_params)*sym_params;
    fprintf(1,'Regrouped parameter %s <= ', char(sym_params(ind)));
    disp(sym_result)
end
```
In this instance, `independent_param_inds` stores the indices of all identifiable inertial parameters.

Note that the dynamics vector is stored in different orders between [spatial_v2_extended](https://github.com/ROAM-Lab-ND/spatial_v2_extended/) and [pinocchio](https://github.com/stack-of-tasks/pinocchio/)!

Order in [spatial_v2_extended](https://github.com/ROAM-Lab-ND/spatial_v2_extended/) can be found in [this line](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/inertia/inertiaMatToVec.m#L8).
```
[m hx hy hz Ixx Iyy Izz Iyz Ixz Ixy]
```
Order in [pinocchio](https://github.com/stack-of-tasks/pinocchio/) can be found in [this line](https://github.com/stack-of-tasks/pinocchio/blob/master/include/pinocchio/spatial/inertia.hpp#L536).
```
[m hx hy hz Ixx Ixy Iyy Ixz Iyz Izz]
```
To use the results in [RPNA_URDF_Example.m](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/identifiability/RPNA_URDF_Example.m) in **RAPTOR**, which uses pinocchio, we would suggest the following code to change the order of the `independent_param_inds`:

```MATLAB
inertia_reorder = [1, 2, 3, 4, 5, 7, 10, 9, 8, 6];
for i = 1:length(independent_param_inds)
    suborder = mod(independent_param_inds(i) - 1, 10);
    link_index = floor((independent_param_inds(i) - 1) / 10);
    if mod(independent_param_inds(i) - 1, 10) <= 3
        fprintf('%d, ', independent_param_inds(i) - 1);
    else
        fprintf('%d, ', 10 * link_index + inertia_reorder(suborder + 1) - 1);
    end
end
fprintf('\n');
```
This gives the correct indices in pinocchio order and in C++ format (0-indexed).
Users can copy the printed array and paste it in `KinovaRegressorExample.cpp`.

### Payload Exciting Trajectory Generation
`KinovaRegressorExamplePayload.cpp` implements a trajectory optimization problem that minimizes the condition number of the last 10 columns of the regressor matrix, which correspond to the inertial parameters of the robot's last link. 
This matrix is usually full-rank.
This is particularly useful for identifying a robotic arm that has gone through the system identification process, but then modified with a new gripper or other unknown payloads at the end effector.

## References

[1] B. Armstrong, "On finding 'exciting' trajectories for identification experiments involving systems with non-linear dynamics," Proceedings. 1987 IEEE International Conference on Robotics and Automation, Raleigh, NC, USA, 1987, pp. 1131-1139, doi: 10.1109/ROBOT.1987.

[2] M. Gautier and W. Khalil, "Exciting trajectories for the identification of base inertial parameters of robots," [1991] Proceedings of the 30th IEEE Conference on Decision and Control, Brighton, UK, 1991, pp. 494-499 vol.1, doi: 10.1109/CDC.1991.261353.

[3] K. Ayusawa, A. Rioux, E. Yoshida, G. Venture and M. Gautier, "Generating persistently exciting trajectory based on condition number optimization," 2017 IEEE International Conference on Robotics and Automation (ICRA), Singapore, 2017, pp. 6518-6524, doi: 10.1109/ICRA.2017.7989770.

[4] W. Rackl, R. Lampariello and G. Hirzinger, "Robot excitation trajectories for dynamic parameter estimation using optimized B-splines," 2012 IEEE International Conference on Robotics and Automation, Saint Paul, MN, USA, 2012, pp. 2042-2047, doi: 10.1109/ICRA.2012.6225279.

[5] G. Venture, K. Ayusawa and Y. Nakamura, "A numerical method for choosing motions with optimal excitation properties for identification of biped dynamics - An application to human," 2009 IEEE International Conference on Robotics and Automation, Kobe, Japan, 2009, pp. 1226-1231, doi: 10.1109/ROBOT.2009.5152264.

[6] S. Farsoni, F. Ferraguti, and M. Bonfè, “Safety-oriented robot payload identification using collision-free path planning and decoupling motions,” Robotics and Computer-Integrated Manufacturing, vol. 59, pp. 189–200, 2019.

[7] P. M. Wensing, G. Niemeyer, and J. J. E. Slotine, “A geometric characterization of observability in inertial parameter identification,” The International Journal of Robotics Research, vol. 43, no. 14, pp. 2274–2302, 2024.