# Installation of required packages

First clone the repository
```
git clone https://github.com/roahmlab/RAPTOR.git
```

## Install Through Docker (Recommended)

You shouldcomplete HSL steps BEFORE you build the docker image otherwise you will have error.
<!-- We choose [HSL](https://www.hsl.rl.ac.uk/) to solve large linear systems in the nonlinear optimization problem.
Check this [github repository](https://github.com/coin-or-tools/ThirdParty-HSL) out and follow the instructions there. -->
<!-- To be more specific, you need to  -->
### HSL
We have selected [HSL](https://www.hsl.rl.ac.uk/) to solve large linear systems in the nonlinear optimization problem. Please follow the instructions below to complete the installation:
1. Compared to the official HSL, we used [ThirdParty-HSL](https://github.com/coin-or-tools/ThirdParty-HSL), which is specifically tailored for COIN-OR projects, particularly Ipopt, offering easier integration and installation.
    ```
    git clone https://github.com/coin-or-tools/ThirdParty-HSL.git. 
    ```
2. Download a tarball containing the Coin-HSL source code from the official [website](https://licences.stfc.ac.uk/product/coin-hsl). The academic license is free but it could take 1 or 2 days to process the order.
3. Unpack the Coin-HSL source code and place the folder inside [ThirdParty-HSL](https://github.com/coin-or-tools/ThirdParty-HSL), which serves as a wrapper to simplify the compilation and integration of HSL.
4. Rename the ThirdParty-HSL as `HSL`,then compress it into a `HSL.zip` file. Note: This zip file will later be 'unziped' and built inside the docker container. <!-- Not sure we should keep one uncompressed copy at origin location -->


### Docker
We strongly recommend using Docker. We have provided a Dockerfile that will automatically install all the required packages. If you don't have Docker installed, you can find the installation instructions [here](https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository).

1. The provided Dockerfile is at `docker/Dockerfile`. 
2. Move the HSL.zip file in `docker/` so that dockerfile can find it.
3. In Visual Studio Code, simply click "ctrl+shift+P" and search "Dev Containers: Rebuild and Reopen Container", it will build the environment automatically for you from [docker/Dockerfile](../docker/Dockerfile) so that you don't need to follow the steps below.

## Install Required Packages Independently
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
sudo apt-get install libboost-all-dev
```

### urdfdom
```
sudo apt-get install liburdfdom-dev
```

### Eigen 3.4
In Ubuntu 22.04, the Eigen library version is 3.4.0 by default. Simply do
```
sudo apt install libeigen3-dev
```
In Ubuntu 20.04, you need to go to the official [website](https://eigen.tuxfamily.org/index.php?title=3.4) and manually install it.

### pinocchio
Please refer to the offical website [link](https://stack-of-tasks.github.io/pinocchio/download.html).

<!-- ### Qhull (Not used for now)
Please refer to this [link](http://www.qhull.org/download/)
Build from source in the downloaded folder so that the libraries are installed in /usr/local/ -->

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