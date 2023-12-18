# Installation of required packages

All of the following instructions are assuming that the packages are installed to root directories,
like /usr/local or /usr.
You might need to edit your ~/.basrc so that the environment variable PATH includes /usr/local/bin and
LD_LIBRARY_PATH includes /usr/local/lib.
Be careful if you are using a shared server.

### GSL
```
sudo apt-get install libgsl-dev
```

### Boost
```
sudo apt-get install libboost-dev
```

### urdfdom
```
sudo apt-get install liburdfdom-dev
```

### PkgConfig
```
sudo apt-get install -y pkg-config
```

### Eigen 3.4
In Ubuntu 22.04, the Eigen library version is 3.4.0 by default. Simply do
```
sudo apt install libeigen3-dev
```
In Ubuntu 20.04, you need to go to the official [website](https://eigen.tuxfamily.org/index.php?title=3.4) and manually install it.

### pinocchio
Please refer to this [link](https://stack-of-tasks.github.io/pinocchio/download.html).
Go to "Build from source" page. 
Follow the instructions there.
Note that if you encounter any problem regarding python binding (like, eigenpy), go to the CMakeLists.txt in the main folder, change flag BUILD_PYTHON_INTERFACE to OFF, and redo all the procedures.
Specifically, run the following command
```
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_PYTHON_INTERFACE=OFF
```

### Ipopt
Please refer to this [link](https://coin-or.github.io/Ipopt/INSTALL.html).

Note that in section "Download, build, and install dependencies", you actually only to install one of the external libraries.
We choose [HSL](https://www.hsl.rl.ac.uk/) to solve large linear systems in a nonlinear optimization problem.

Check this [github repository](https://github.com/coin-or-tools/ThirdParty-HSL) out and follow the instructions there.
To be more specific, you need to get HSL library from its official [website](https://www.hsl.rl.ac.uk/), 
and then put the folder inside [ThirdParty-HSL](https://github.com/coin-or-tools/ThirdParty-HSL).
After compilation, you should be able to find the library "libcoinhsl.so" installed in /usr/local/lib.
If not, please manually move them to /usr/local/lib so that it is easier for cmake to find them automatically.

I have had one issue before where ipopt tries to find "libhsl.so" instead of "libcoinhsl.so". 
So "libhsl.so" is the actual name of this library while "libcoinhsl.so" is the name when compiled in the [ThirdParty-HSL](https://github.com/coin-or-tools/ThirdParty-HSL) wrapper.
You might need to create a symbolic link manually 
```
sudo ln -s libcoinhsl.so libhsl.so
```
There are a lot of moving pieces here so something might not work or my README is still not enough.
Please contact me (jimzhang@umich.edu) if you encounter any problems installing ipopt.