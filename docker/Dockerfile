# Use an official Ubuntu runtime as a parent image
FROM nvidia/cuda:12.1.1-devel-ubuntu22.04

# Set the working directory in the container
WORKDIR /app

SHELL ["/bin/bash", "-c"]

# Copy the HSL source code folder into the Docker image
COPY docker/HSL.zip /app
# COPY docker/hsl_ma57.zip /app
# COPY docker/hsl_ma86.zip /app

# # Copy the pardiso libraries into the Docker image
# COPY docker/pardiso/include/pardiso.h /usr/local/include/pardiso.h
# COPY docker/pardiso/lib/libpardiso.so /usr/local/lib/libpardiso.so
# COPY docker/pardiso/lib/libpardiso_no_deps.so /usr/local/lib/libpardiso_no_deps.so
# COPY docker/pardiso/panua.lic /app/panua.lic
# ENV  PARDISO_LIC_PATH /app

# Install necessary packages
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    gdb \
    valgrind \
    libboost-all-dev \
    libopenblas-dev \
    libeigen3-dev \
    libgsl-dev \
    git \
    unzip \
    patch \
    wget \
    pkg-config \
    liblapack-dev \
    libmetis-dev \
    python3.10 \
    python3-pip \
    iputils-ping \
    net-tools \
    netcat \
    screen \
    lsb-release \
    gnupg2 \
    curl \
    libyaml-cpp-dev \
    xvfb

#     # install python packages
# RUN pip install --upgrade pip && \
#     pip install scipy matplotlib pyyaml torch pybullet==3.2.5 mujoco==2.3.7 glfw==2.6.2 pyopengl==3.1.7 numpy==1.25.2 nanobind pygccxml pyplusplus
RUN pip install --upgrade pip && \
    pip install scipy matplotlib pyyaml pybullet==3.2.5 mujoco==2.3.7 glfw==2.6.2 pyopengl==3.1.7 numpy==1.25.2 nanobind pygccxml pyplusplus

    # Install pinochio dependencies
RUN --mount=type=cache,target=/var/cache/apt \
    --mount=type=bind,source=docker/scripts/install-pinocchio.sh,target=/tmp/install-pinocchio.sh \
    bash /tmp/install-pinocchio.sh

    # Install the HSL solver
RUN --mount=type=cache,target=/var/cache/apt \
    --mount=type=bind,source=docker/scripts/install-hsl.sh,target=/tmp/install-hsl.sh \
    bash /tmp/install-hsl.sh

    # Install the Ipopt solver
RUN --mount=type=cache,target=/var/cache/apt \
    --mount=type=bind,source=docker/scripts/install-ipopt.sh,target=/tmp/install-ipopt.sh \
    bash /tmp/install-ipopt.sh

# Install the ompl
RUN --mount=type=cache,target=/var/cache/apt \
    --mount=type=bind,source=docker/scripts/install-ompl.sh,target=/tmp/install-ompl.sh \
    bash /tmp/install-ompl.sh

# Download and extract MOSEK
# COPY docker/scripts/mosektoolslinux64x86.tar.bz2 /tmp/mosek.tar.bz2
# RUN mkdir -p /opt/mosek && \
#     tar -xjf /tmp/mosek.tar.bz2 -C /opt/mosek --strip-components=1 && \
#     rm /tmp/mosek.tar.bz2

# COPY docker/scripts/mosek.lic /opt/mosek/mosek.lic

# ENV PATH=/opt/mosek/10.2/tools/platform/linux64x86/bin:$PATH
# ENV MOSEKLM_LICENSE_FILE=/opt/mosek/mosek.lic

# # MOSEK env
# ENV MOSEK_HOME=/opt/mosek/10.2
# ENV PATH=$MOSEK_HOME/tools/platform/linux64x86/bin:$PATH
# ENV LD_LIBRARY_PATH=$MOSEK_HOME/tools/platform/linux64x86/bin:$LD_LIBRARY_PATH
# ENV CMAKE_PREFIX_PATH=$MOSEK_HOME/tools/platform/linux64x86/h:$CMAKE_PREFIX_PATH 

# # # Compile MOSEK Fusion C++ API
# WORKDIR $MOSEK_HOME/tools/platform/linux64x86/src/fusion_cxx
# RUN make install   


ENV PATH /opt/openrobots/bin:/usr/local/bin:/usr/bin:$PATH
ENV PKG_CONFIG_PATH /opt/openrobots/lib/pkgconfig
ENV LD_LIBRARY_PATH /usr/lib:/usr/lib/x86_64-linux-gnu:/opt/openrobots/lib:/usr/local/lib:/usr/local/share:$LD_LIBRARY_PATH
ENV PYTHONPATH /opt/openrobots/lib/python3.10/site-packages
ENV CMAKE_PREFIX_PATH /usr:/opt/openrobots:/usr/local:/usr/local/include/ompl-1.6:$CMAKE_PREFIX_PATH 

# Set the entry point to VSCode's default command
CMD ["/bin/bash"]

USER $USERNAME