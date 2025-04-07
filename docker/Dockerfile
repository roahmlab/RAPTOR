# Use an official Ubuntu runtime as a parent image
FROM ubuntu:22.04

# Set the working directory in the container
WORKDIR /app

SHELL ["/bin/bash", "-c"]

# Copy the HSL source code folder into the Docker image
COPY docker/HSL.zip /app

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
RUN pip install --upgrade pip && \
    pip install scipy matplotlib pyyaml pybullet==3.2.5 glfw==2.6.2 pyopengl==3.1.7 numpy==1.25.2 nanobind pygccxml pyplusplus

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

ENV PATH /opt/openrobots/bin:/usr/local/bin:/usr/bin:$PATH
ENV PKG_CONFIG_PATH /opt/openrobots/lib/pkgconfig
ENV LD_LIBRARY_PATH /usr/lib:/usr/lib/x86_64-linux-gnu:/opt/openrobots/lib:/usr/local/lib:/usr/local/share:$LD_LIBRARY_PATH
ENV PYTHONPATH /opt/openrobots/lib/python3.10/site-packages
ENV CMAKE_PREFIX_PATH /usr:/opt/openrobots:/usr/local:/usr/local/include/ompl-1.7:$CMAKE_PREFIX_PATH 

# Set the entry point to VSCode's default command
CMD ["/bin/bash"]

USER $USERNAME