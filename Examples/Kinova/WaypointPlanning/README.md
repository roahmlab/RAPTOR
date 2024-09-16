# Waypoints Planning

This example is not involved with optimization.
We incorporate the collision checking functions defined in class [KinovaCustomizedConstraints](KinovaCustomizedConstraints.h) into [ompl](https://ompl.kavrakilab.org/) to plan for a collision-free path (not a trajectory) using RRT-like methods.