# Installation of required packages

First clone the repository
```
git clone https://github.com/roahmlab/RAPTOR.git
```

## Install Through Docker (Strongly Recommended)

We strongly recommend using Docker. We have provided a Dockerfile that will automatically install all the required packages. If you don't have Docker installed, you can find the installation instructions [here](https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository).

### HSL
You should complete HSL steps BEFORE you build the docker image otherwise you will have error.

We have selected [HSL](https://www.hsl.rl.ac.uk/) to solve large linear systems in the nonlinear optimization problem. 
Please follow the instructions below to complete the installation.

Besides official HSL code, we used [ThirdParty-HSL](https://github.com/coin-or-tools/ThirdParty-HSL), which is specifically tailored for COIN-OR projects, particularly Ipopt, offering easier integration and installation.
    ```
    git clone https://github.com/coin-or-tools/ThirdParty-HSL.git. 
    ```
1. Download the tarball containing the Coin-HSL source code from its official [website](https://licences.stfc.ac.uk/product/coin-hsl). You will need to apply a license for it. The academic license is free but it could take 1 or 2 days to process the order.
2. Download code from [ThirdParty-HSL](https://github.com/coin-or-tools/ThirdParty-HSL).
3. Unpack the Coin-HSL source code and rename the folder as `coinhsl`.
4. Place the `coinhsl` folder inside [ThirdParty-HSL](https://github.com/coin-or-tools/ThirdParty-HSL) folder, which serves as a wrapper to simplify the compilation and integration of HSL.
5. Rename the [ThirdParty-HSL](https://github.com/coin-or-tools/ThirdParty-HSL) folder as `HSL`,then compress it into a `HSL.zip` file. Note: This zip file will later be 'unziped' and built inside the docker container.

### Docker
The provided Dockerfile is at [docker/Dockerfile](../docker/Dockerfile). 
Follow the following intructions to build and enter the docker image.

1. Move the `HSL.zip` file in `docker/` so that dockerfile can find it.
2. In Visual Studio Code, simply click `ctrl+shift+P` and search "Dev Containers: Rebuild and Reopen Container", it will build the environment automatically for you from [docker/Dockerfile](../docker/Dockerfile) so that you don't need to follow the steps below.

## Install Required Packages Independently
All of the following instructions are assuming that the packages are installed to root directories,
like /usr/local or /usr.
You might need to edit your ~/.basrc so that the environment variable PATH includes /usr/local/bin and
LD_LIBRARY_PATH includes /usr/local/lib.
Be careful if you are using a shared server.

### GSL
```shell
sudo apt-get install libgsl-dev
```

### Boost
```shell
sudo apt-get install libboost-all-dev
```

### yaml-cpp
```shell
sudo apt-get install libyaml-cpp-dev
```

<!-- ### urdfdom
```shell
sudo apt-get install liburdfdom-dev
``` -->

### Eigen 3.4
In Ubuntu 22.04, the Eigen library version is 3.4.0 by default. Simply do
```shell
sudo apt install libeigen3-dev
```
In Ubuntu 20.04, the Eigen library version is 3.3.7 by default.
You will have to go to the official [website](https://eigen.tuxfamily.org/index.php?title=3.4) and manually install it.

### pinocchio
Please refer to the offical website [link](https://stack-of-tasks.github.io/pinocchio/download.html).

One small issue could be: `pinocchio` would recommend you install through `robotpkg` on Linux systems, while some people mainly works with `pinocchio` in Python so they install `pinocchio` through `conda`.
`robotpkg` installs `pinocchio` in `/opt/openrobots/` while `conda` installs `pinocchio` in `~/miniconda3/`.
Try to avoid installing two `pinocchio` on your computer.
If you have to, be sure to specify the paths carefully when you try to link your program to `pinocchio`.

### Ipopt
Please refer to this [link](https://coin-or.github.io/Ipopt/INSTALL.html).

Note that in section "Download, build, and install dependencies", you actually only to install one of the external libraries.
We choose [HSL](https://www.hsl.rl.ac.uk/) to solve large linear systems in a nonlinear optimization problem.

Check this [github repository](https://github.com/coin-or-tools/ThirdParty-HSL) out and follow the instructions there.
To be more specific, you need to get HSL library from its official [website](https://www.hsl.rl.ac.uk/), 
and then put the folder inside [ThirdParty-HSL](https://github.com/coin-or-tools/ThirdParty-HSL).
After compilation, you should be able to find the library `libcoinhsl.so` installed in /usr/local/lib.
If not, please manually move them to `/usr/local/lib` so that it is easier for cmake to find them automatically.

I have had one issue before where ipopt tries to find `libhsl.so` instead of `libcoinhsl.so`. 
So `libhsl.so` is the actual name of this library while `libcoinhsl.so` is the name when compiled in the [ThirdParty-HSL](https://github.com/coin-or-tools/ThirdParty-HSL) wrapper.
You might need to create a symbolic link manually 
```shell
sudo ln -s libcoinhsl.so libhsl.so
```
There are a lot of moving pieces here so something might not work or my README is still not enough.
Please raise a github [issue](https://github.com/roahmlab/RAPTOR/issues) if you encounter any problems installing ipopt.