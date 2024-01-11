# Minimum required CMake version
cmake_minimum_required(VERSION 3.12)

# Project name
project(IDTO
    DESCRIPTION "Inverse Dynamics Trajectory Optimization")

# if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release)
# endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp")

message(${CMAKE_BUILD_TYPE})

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)

# Find Eigen library
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

# Find GSL library
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

# Find Qhull library
find_package(Qhull REQUIRED)
include_directories(${QHULL_INCLUDE_DIRS})

# Add pinocchio configurations
include(cmake/pinocchio_config.cmake)

# Find urdf library
#find_package(Urdfdom REQUIRED COMPONENTS
#    urdf
#)
#include_directories(${URDFDOM_INCLUDE_DIRS})

# Find Boost library
find_package(Boost REQUIRED COMPONENTS
    system
    filesystem
)

# Find Ipopt library
set(IPOPT_LIBRARY_DIR /usr/local/lib)
set(IPOPT_INCLUDE_DIR /usr/local/include/coin-or)
include_directories(${IPOPT_INCLUDE_DIR})
link_directories(${IPOPT_LIBRARY_DIR})

# Find HSL library
set(HSL_LIBRARY_DIR /usr/local/lib)
set(HSL_INCLUDE_DIR /usr/local/include)
include_directories(${HSL_INCLUDE_DIR})
link_directories(${HSL_LIBRARY_DIR})

find_package(PkgConfig REQUIRED)
pkg_check_modules(PKG_PIN_CONFIG REQUIRED 
    pinocchio)
include_directories(${PKG_PIN_CONFIG_INCLUDE_DIRS})

## Trajectories library
add_library(trajlib SHARED
    Trajectories/src/FourierCurves.cpp 
    Trajectories/src/FixedFrequencyFourierCurves.cpp
    Trajectories/src/BezierCurves.cpp 
    Trajectories/src/Polynomials.cpp 
    Trajectories/src/Plain.cpp 
    Trajectories/src/Trajectories.cpp)

target_include_directories(trajlib PUBLIC
    Trajectories/include
    Utils/include)

## Inverse dynamics library
add_library(IDlib SHARED
    InverseDynamics/src/InverseDynamics.cpp
    InverseDynamics/src/DynamicsConstraints.cpp
    InverseDynamics/src/ConstrainedInverseDynamics.cpp
    InverseDynamics/src/Transform.cpp
    InverseDynamics/src/ForwardKinematics.cpp)

target_include_directories(IDlib PUBLIC
    InverseDynamics/include
    Trajectories/include
    Utils/include)

target_compile_options(IDlib PUBLIC
    ${PINOCCHIO_FLAGS})

## Constraints library
add_library(Conslib SHARED
    Constraints/src/JointLimits.cpp
    Constraints/src/ConstrainedJointLimits.cpp
    Constraints/src/TorqueLimits.cpp
    Constraints/src/SurfaceContactConstraints.cpp
    Constraints/src/KinematicsConstraints.cpp
    Constraints/src/EndEffectorConstraints.cpp
    Constraints/src/ZonotopeCollisionAvoidance.cpp)

target_include_directories(Conslib PUBLIC
    InverseDynamics/include
    Trajectories/include
    Constraints/include
    Utils/include
    /usr/local/include/libqhullcpp)

target_compile_options(Conslib PUBLIC
    ${PINOCCHIO_FLAGS})    

## Constraints library
add_library(Optlib SHARED
    Optimization/src/Optimizer.cpp)

target_include_directories(Optlib PUBLIC
    Optimization/include
    Constraints/include
    Utils/include)

## Digit library
add_library(Digitlib SHARED
    Examples/Digit/src/DigitDynamicsConstraints.cpp
    Examples/Digit/src/DigitConstrainedInverseDynamics.cpp
    Examples/Digit/src/DigitSingleStepOptimizer.cpp
    Examples/Digit/src/DigitCustomizedConstraints.cpp
    Examples/Digit/src/DigitSingleStepPeriodicityConstraints.cpp)

target_link_libraries(Digitlib PRIVATE
    trajlib
    IDlib
    Conslib
    Optlib
    ipopt
    coinhsl
    ${EIGEN3_LIBRARIES} 
    ${GSL_LIBRARIES}
    ${BOOST_LIBRARIES} 
    ${URDFDOM_LIBRARIES}
    ${PKG_PIN_CONFIG_LIBRARIES})

target_include_directories(Digitlib PUBLIC
    Examples/Digit/include
    InverseDynamics/include
    Trajectories/include
    Constraints/include
    Optimization/include
    Utils/include)

target_compile_options(Digitlib PUBLIC
    ${PINOCCHIO_FLAGS})

## Digit-modified library
add_library(DigitModifiedlib SHARED
    Examples/Digit-modified/src/DigitModifiedDynamicsConstraints.cpp
    Examples/Digit-modified/src/DigitModifiedConstrainedInverseDynamics.cpp
    Examples/Digit-modified/src/DigitModifiedSingleStepOptimizer.cpp
    Examples/Digit-modified/src/DigitModifiedCustomizedConstraints.cpp
    Examples/Digit-modified/src/DigitModifiedSingleStepPeriodicityConstraints.cpp)

target_link_libraries(DigitModifiedlib PRIVATE
    trajlib
    IDlib
    Conslib
    Optlib
    ipopt
    coinhsl
    ${EIGEN3_LIBRARIES} 
    ${BOOST_LIBRARIES} 
    ${URDFDOM_LIBRARIES}
    ${PKG_PIN_CONFIG_LIBRARIES})

target_include_directories(DigitModifiedlib PUBLIC
    Examples/Digit-modified/include
    InverseDynamics/include
    Trajectories/include
    Constraints/include
    Optimization/include
    Utils/include)

target_compile_options(DigitModifiedlib PUBLIC
    ${PINOCCHIO_FLAGS})

# Kinova library
add_library(Kinovalib SHARED
    Examples/Kinova/src/KinovaOptimizer.cpp
    Examples/Kinova/src/KinovaCustomizedConstraints.cpp)

target_link_libraries(Kinovalib PRIVATE
    trajlib
    IDlib
    Conslib
    Optlib
    ipopt
    coinhsl
    ${EIGEN3_LIBRARIES} 
    ${BOOST_LIBRARIES} 
    ${URDFDOM_LIBRARIES}
    ${PKG_PIN_CONFIG_LIBRARIES}
    ${QHULL_LIBRARIES})

target_include_directories(Kinovalib PUBLIC
    Examples/Kinova/include
    InverseDynamics/include
    Trajectories/include
    Constraints/include
    Optimization/include
    Utils/include)

target_compile_options(Kinovalib PUBLIC
    ${PINOCCHIO_FLAGS})

## Example scripts
add_executable(Digit_example 
    Examples/Digit/DigitSingleStep.cpp)

target_include_directories(Digit_example PRIVATE
    Examples/Digit/include
    InverseDynamics/include
    Trajectories/include
    Constraints/include
    Optimization/include
    Utils/include)

target_link_libraries(Digit_example PRIVATE
    trajlib 
    IDlib
    Conslib
    Optlib
    Digitlib
    ipopt
    coinhsl
    ${EIGEN3_LIBRARIES} 
    ${GSL_LIBRARIES}
    ${BOOST_LIBRARIES} 
    ${URDFDOM_LIBRARIES}
    ${PKG_PIN_CONFIG_LIBRARIES})

target_compile_options(Digit_example PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(DigitModified_example 
    Examples/Digit-modified/DigitModifiedSingleStep.cpp)

target_include_directories(DigitModified_example PRIVATE
    Examples/Digit-modified/include
    InverseDynamics/include
    Trajectories/include
    Constraints/include
    Optimization/include
    Utils/include)

target_link_libraries(DigitModified_example PRIVATE
    trajlib 
    IDlib
    Conslib
    Optlib
    DigitModifiedlib
    ipopt
    coinhsl
    ${EIGEN3_LIBRARIES} 
    ${BOOST_LIBRARIES} 
    ${URDFDOM_LIBRARIES}
    ${PKG_PIN_CONFIG_LIBRARIES})

target_compile_options(DigitModified_example PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(Kinova_example 
    Examples/Kinova/KinovaExample.cpp)

target_include_directories(Kinova_example PRIVATE
    Examples/Kinova/include
    InverseDynamics/include
    Trajectories/include
    Constraints/include
    Optimization/include
    Utils/include)

target_link_libraries(Kinova_example PRIVATE
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${EIGEN3_LIBRARIES} 
    ${BOOST_LIBRARIES} 
    ${URDFDOM_LIBRARIES}
    ${PKG_PIN_CONFIG_LIBRARIES})

target_compile_options(Kinova_example PUBLIC
    ${PINOCCHIO_FLAGS})

## Install
install(TARGETS trajlib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(TARGETS IDlib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(TARGETS Conslib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(TARGETS Optlib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(TARGETS Digitlib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(TARGETS Kinovalib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(TARGETS Digit_example
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(TARGETS DigitModified_example
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(TARGETS Kinova_example
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
