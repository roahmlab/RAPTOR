# Minimum required CMake version
cmake_minimum_required(VERSION 3.12)

# Project name
project(IDTO
    DESCRIPTION "Inverse Dynamics Trajectory Optimization")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)

# Find Eigen library
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

# Find GSL library
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

# Add pinocchio configurations
include(cmake/pinocchio_config.cmake)

# Find urdf library
find_package(Urdfdom REQUIRED COMPONENTS
    urdf
)
include_directories(${URDFDOM_INCLUDE_DIRS})

# Find Boost library
find_package(Boost REQUIRED COMPONENTS
    system
    filesystem
)

# Find Ipopt library
set(IPOPT_LIBRARIES /usr/local/lib)
set(IPOPT_INCLUDE_DIR /usr/local/include/coin-or)
include_directories(${IPOPT_INCLUDE_DIR})

find_package(PkgConfig REQUIRED)
pkg_check_modules(PKG_PIN_CONFIG REQUIRED 
    pinocchio)
include_directories(${PKG_PIN_CONFIG_INCLUDE_DIRS})

## Trajectories library
add_library(trajlib SHARED
    Trajectories/src/FourierCurves.cpp 
    Trajectories/src/Trajectories.cpp)

target_include_directories(trajlib PUBLIC
    Trajectories/include)

## Inverse dynamics library
add_library(IDlib SHARED
    InverseDynamics/src/InverseDynamics.cpp
    InverseDynamics/src/DynamicsConstraints.cpp
    InverseDynamics/src/ConstrainedInverseDynamics.cpp
    InverseDynamics/src/Transform.cpp
    InverseDynamics/src/ForwardKinematics.cpp)

target_include_directories(IDlib PUBLIC
    InverseDynamics/include)

target_compile_options(IDlib PUBLIC
    ${PINOCCHIO_FLAGS})

## Digit library
add_library(Digitlib SHARED
    Examples/Digit/src/DigitDynamicsConstraints.cpp
    Examples/Digit/src/DigitConstrainedInverseDynamics.cpp
    Examples/Digit/src/DigitSingleStepOptimizer.cpp)

target_link_libraries(Digitlib PRIVATE
    trajlib
    IDlib
    ${GSL_LIBRARIES})

target_include_directories(Digitlib PUBLIC
    Examples/Digit/include
    InverseDynamics/include)

target_compile_options(Digitlib PUBLIC
    ${PINOCCHIO_FLAGS})

## Example scripts
add_executable(IDTO_example 
    Examples/Digit/DigitSingleStep.cpp)

target_include_directories(IDTO_example PRIVATE
    Trajectories/include)

target_link_libraries(IDTO_example PRIVATE
    trajlib 
    IDlib
    Digitlib
    ${EIGEN3_LIBRARIES} 
    ${GSL_LIBRARIES}
    ${BOOST_LIBRARIES} 
    ${IPOPT_LIBRARIES}
    ${URDFDOM_LIBRARIES}
    ${PKG_PIN_CONFIG_LIBRARIES})

target_compile_options(IDTO_example PUBLIC
    ${PINOCCHIO_FLAGS})

## Install
install(TARGETS trajlib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(TARGETS IDlib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(TARGETS IDTO_example
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
