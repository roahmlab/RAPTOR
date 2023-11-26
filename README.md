# IDTO
Trajectory Optimization based on Inverse Dynamics

This branch attempts to implement the sparse gradient for Ipopt. 
But it is unfinished right now and is not verified to be useful (significantly faster than dense gradient).

## Requirement
- Ubuntu >= 20.04
- [Eigen 3.4](https://eigen.tuxfamily.org/index.php?title=3.4)
- [pinocchio](https://stack-of-tasks.github.io/pinocchio/download.html)
- [Ipopt](https://coin-or.github.io/Ipopt/INSTALL.html)
- Boost (libboost-dev)
- GSL (libgsl-dev)

## Getting Started
Run the following commands for Digit Example
```
mkdir build
cd build
cmake ..
make
./IDTO_example
```

## Credits
Bohao Zhang (jimzhang@umich.edu)

roahmlab, Umich
