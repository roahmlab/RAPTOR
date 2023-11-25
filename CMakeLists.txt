# Minimum required CMake version
cmake_minimum_required(VERSION 3.12)

# Project name
project(IDTO
    DESCRIPTION "Inverse Dynamics Trajectory Optimization")

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
    Examples/Digit/src/DigitConstrainedInverseDynamics.cpp)

target_link_libraries(Digitlib PRIVATE
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
    ${GSL_LIBRARIES})

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
