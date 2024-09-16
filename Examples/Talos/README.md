# Gait Optimization For Talos

This folder contains the gait optimization implementation for Talos, that is used for comparison between RAPTOR and [aligator](https://github.com/Simple-Robotics/aligator/).

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