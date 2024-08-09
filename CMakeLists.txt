# Minimum required CMake version
cmake_minimum_required(VERSION 3.18)

# Project name
project(RAPTOR
    DESCRIPTION "Rapid and Robust Trajectory Optimization for Robots")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp")

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Number of threads that openmp will use
if(DEFINED NUM_THREADS)
  add_definitions(-DNUM_THREADS=${NUM_THREADS})
  message(STATUS "Number of threads used in openmp: ${NUM_THREADS}")
else()
    # Get the maximum number of threads using nproc
    execute_process(
        COMMAND nproc
        OUTPUT_VARIABLE MAX_THREADS
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    add_definitions(-DNUM_THREADS=${MAX_THREADS})
    message(STATUS "Maximize number of threads used in openmp: ${MAX_THREADS}")
endif()

# Find Python binding library
find_package(Python 3.10
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# Find Boost library
find_package(Boost REQUIRED COMPONENTS 
    system
    filesystem
    serialization
)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

# Find GSL library
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})

# Find Eigen library
find_package(Eigen3 3.3.7 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(pinocchio REQUIRED)
include_directories(${PINOCCHIO_INCLUDE_DIRS})
link_directories(${PINOCCHIO_LIBRARY_DIRS})

set(PINOCCHIO_FLAGS
  ${pinocchio_CFLAGS_OTHER}
  -Wno-ignored-attributes
  -Wno-invalid-partial-specialization    # to silence warning with unsupported Eigen Tensor
  -DPINOCCHIO_URDFDOM_TYPEDEF_SHARED_PTR
  -DPINOCCHIO_URDFDOM_USE_STD_SHARED_PTR
)

# Find YAML library
find_package(yaml-cpp REQUIRED)
include_directories(${YAML_CPP_INCLUDE_DIR})
link_directories(${YAML_CPP_LIBRARY_DIR})

# Find Ipopt library
set(IPOPT_INCLUDE_DIR /usr/local/include/coin-or)
set(IPOPT_LIBRARY_DIR /usr/local/lib)
include_directories(${IPOPT_INCLUDE_DIR})
link_directories(${IPOPT_LIBRARY_DIR})

# Find HSL library
set(HSL_INCLUDE_DIR /usr/local/include/coin/ThirdParty)
set(HSL_LIBRARY_DIR /usr/local/lib)
include_directories(${HSL_INCLUDE_DIR})
link_directories(${HSL_LIBRARY_DIR})

# Find pardiso library
set(PARDISO_INCLUDE_DIR /usr/local/include)
set(PARDISO_LIBRARY_DIR /usr/local/lib)
include_directories(${PARDISO_INCLUDE_DIR})
link_directories(${PARDISO_LIBRARY_DIR})

include_directories(Utils/include)

## Trajectories library
add_library(trajlib SHARED
    Trajectories/src/Trajectories.cpp
    Trajectories/src/FourierCurves.cpp 
    Trajectories/src/FixedFrequencyFourierCurves.cpp
    Trajectories/src/BezierCurves.cpp 
    Trajectories/src/ArmourBezierCurves.cpp
    Trajectories/src/Polynomials.cpp 
    Trajectories/src/Plain.cpp
    Trajectories/src/TrajectoryGroup.cpp)

target_include_directories(trajlib PUBLIC
    Trajectories/include)

include_directories(Trajectories/include)

## Inverse dynamics library
add_library(IDlib SHARED
    KinematicsDynamics/src/InverseDynamics.cpp
    KinematicsDynamics/src/CustomizedInverseDynamics.cpp
    KinematicsDynamics/src/RegressorInverseDynamics.cpp
    KinematicsDynamics/src/DynamicsConstraints.cpp
    KinematicsDynamics/src/ConstrainedInverseDynamics.cpp
    KinematicsDynamics/src/Transform.cpp
    KinematicsDynamics/src/Spatial.cpp
    KinematicsDynamics/src/ForwardKinematics.cpp)

target_include_directories(IDlib PUBLIC
    KinematicsDynamics/include)

include_directories(KinematicsDynamics/include)

## Constraints library
add_library(Conslib SHARED
    Constraints/src/Constraints.cpp
    Constraints/src/JointLimits.cpp
    Constraints/src/VelocityLimits.cpp
    Constraints/src/ConstrainedJointLimits.cpp
    Constraints/src/TorqueLimits.cpp
    Constraints/src/SurfaceContactConstraints.cpp
    Constraints/src/WaitrContactConstraints.cpp
    Constraints/src/LieSpaceResidual.cpp
    Constraints/src/KinematicsConstraints.cpp
    Constraints/src/BoxCollisionAvoidance.cpp)
    # Constraints/src/ZonotopeCollisionAvoidance.cpp
    # Constraints/src/fclCollisionAvoidance.cpp)

target_include_directories(Conslib PUBLIC
    Constraints/include)
    
include_directories(Constraints/include)

## Constraints library
add_library(Optlib SHARED
    Optimization/src/Optimizer.cpp)

target_include_directories(Optlib PUBLIC
    Optimization/include)

include_directories(Optimization/include)

## Digit library
add_library(Digitlib SHARED
    Examples/Digit/src/DigitDynamicsConstraints.cpp
    Examples/Digit/src/DigitConstrainedInverseDynamics.cpp
    Examples/Digit/src/DigitSingleStepOptimizer.cpp
    Examples/Digit/src/DigitMultipleStepOptimizer.cpp
    Examples/Digit/src/DigitCustomizedConstraints.cpp
    Examples/Digit/src/DigitSingleStepPeriodicityConstraints.cpp
    Examples/Digit/src/DigitMultipleStepPeriodicityConstraints)

target_link_libraries(Digitlib PUBLIC
    trajlib
    IDlib
    Conslib
    Optlib
    ipopt
    coinhsl
    ${GSL_LIBRARIES}
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_include_directories(Digitlib PUBLIC
    Examples/Digit/include)

## Talos library
add_library(Taloslib SHARED
    Examples/Talos/src/TalosDynamicsConstraints.cpp
    Examples/Talos/src/TalosConstrainedInverseDynamics.cpp
    Examples/Talos/src/TalosCustomizedConstraints.cpp
    Examples/Talos/src/TalosSingleStepOptimizer.cpp
    Examples/Talos/src/TalosMultipleStepOptimizer.cpp
    Examples/Talos/src/TalosMultipleStepPeriodicityConstraints)

target_link_libraries(Taloslib PUBLIC
    trajlib
    IDlib
    Conslib
    Optlib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_include_directories(Taloslib PUBLIC
    Examples/Talos/include)

# ## Digit-modified library
# add_library(DigitModifiedlib SHARED
#     Examples/Digit-modified/src/DigitModifiedDynamicsConstraints.cpp
#     Examples/Digit-modified/src/DigitModifiedConstrainedInverseDynamics.cpp
#     Examples/Digit-modified/src/DigitModifiedSingleStepOptimizer.cpp
#     Examples/Digit-modified/src/DigitModifiedCustomizedConstraints.cpp
#     Examples/Digit-modified/src/DigitModifiedSingleStepPeriodicityConstraints.cpp)

# target_link_libraries(DigitModifiedlib PUBLIC
#     trajlib
#     IDlib
#     Conslib
#     Optlib
#     ipopt
#     coinhsl
#     ${BOOST_LIBRARIES} 
#     pinocchio::pinocchio)

# target_include_directories(DigitModifiedlib PUBLIC
#     Examples/Digit-modified/include)

# Kinova libraries
add_library(Kinovalib SHARED
    Examples/Kinova/InverseKinematics/KinovaIKSolver.cpp
    Examples/Kinova/ArmourDiscrete/src/KinovaOptimizer.cpp
    Examples/Kinova/ArmourDiscrete/src/KinovaCustomizedConstraints.cpp
    Examples/Kinova/WaitrDiscrete/src/KinovaWaitrOptimizer.cpp
    Examples/Kinova/SystemIdentification/DataFilter/src/DataFilterOptimizer.cpp
    Examples/Kinova/SystemIdentification/ExcitingTrajectories/src/ConditionNumberOptimizer.cpp)

target_link_libraries(Kinovalib PUBLIC
    trajlib
    IDlib
    Conslib
    Optlib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_include_directories(Kinovalib PUBLIC
    Examples/Kinova
    Examples/Kinova/InverseKinematics
    Examples/Kinova/ArmourDiscrete/include
    Examples/Kinova/WaitrDiscrete/include
    Examples/Kinova/SystemIdentification/DataFilter/include
    Examples/Kinova/SystemIdentification/ExcitingTrajectories/include
    Examples/Kinova/SystemIdentification/IterativeSystemIdentification/include)

## Example scripts
add_executable(Digit_example 
    Examples/Digit/DigitSingleStep.cpp)
    # Examples/Digit/DigitSingleStepSpeedTest.cpp)
    # Examples/Digit/DigitSingleStepRobustnessTest.cpp)
    # Examples/Digit/DigitSingleStepFixedPosition.cpp)

target_include_directories(Digit_example PUBLIC
    Examples/Digit/include
    KinematicsDynamics/include
    Trajectories/include
    Constraints/include
    Optimization/include
    Utils/include)

target_link_libraries(Digit_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Digitlib
    ipopt
    coinhsl
    yaml-cpp
    ${GSL_LIBRARIES}
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(Digit_example PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(Talos_example 
    # Examples/Talos/TalosSingleStep.cpp)
    Examples/Talos/TalosSingleStepFixedPosition.cpp)

target_include_directories(Talos_example PUBLIC
    Examples/Talos/include
    KinematicsDynamics/include
    Trajectories/include
    Constraints/include
    Optimization/include
    Utils/include)

target_link_libraries(Talos_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Taloslib
    ipopt
    coinhsl
    yaml-cpp
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(Talos_example PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(DigitMultipleStep_example 
    Examples/Digit/DigitMultipleStep.cpp)

target_include_directories(DigitMultipleStep_example PRIVATE
    Examples/Digit/include
    KinematicsDynamics/include
    Trajectories/include
    Constraints/include
    Optimization/include
    Utils/include)

target_link_libraries(DigitMultipleStep_example PRIVATE
    trajlib 
    IDlib
    Conslib
    Optlib
    Digitlib
    ipopt
    coinhsl
    yaml-cpp
    ${GSL_LIBRARIES}
    ${BOOST_LIBRARIES}
    pinocchio::pinocchio)

target_compile_options(DigitMultipleStep_example PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(TalosMultipleStep_example 
    Examples/Talos/TalosMultipleStep.cpp)

target_include_directories(TalosMultipleStep_example PRIVATE
    Examples/Talos/include
    KinematicsDynamics/include
    Trajectories/include
    Constraints/include
    Optimization/include
    Utils/include)

target_link_libraries(TalosMultipleStep_example PRIVATE
    trajlib 
    IDlib
    Conslib
    Optlib
    Taloslib
    ipopt
    coinhsl
    yaml-cpp
    ${BOOST_LIBRARIES}
    pinocchio::pinocchio)

target_compile_options(TalosMultipleStep_example PUBLIC
    ${PINOCCHIO_FLAGS})

# add_executable(DigitModified_example 
#     Examples/Digit-modified/DigitModifiedSingleStep.cpp)
#     # Examples/Digit-modified/DigitModifiedSingleStepSpeedTest.cpp)
#     # Examples/Digit-modified/DigitModifiedSingleStepRobustnessTest.cpp)

# target_include_directories(DigitModified_example PUBLIC
#     Examples/Digit-modified/include)

# target_link_libraries(DigitModified_example PUBLIC
#     trajlib 
#     IDlib
#     Conslib
#     Optlib
#     DigitModifiedlib
#     ipopt
#     coinhsl
#     yaml-cpp
#     ${BOOST_LIBRARIES} 
#     pinocchio::pinocchio)

# target_compile_options(DigitModified_example PUBLIC
#     ${PINOCCHIO_FLAGS})

add_executable(KinovaIK_example 
    Examples/Kinova/InverseKinematics/KinovaIKExample.cpp)

target_include_directories(KinovaIK_example PUBLIC
    Examples/Kinova
    Examples/Kinova/InverseKinematics)

target_link_libraries(KinovaIK_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(KinovaIK_example PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(Kinova_example 
    Examples/Kinova/ArmourDiscrete/KinovaExample.cpp)

target_include_directories(Kinova_example PUBLIC
    Examples/Kinova
    Examples/Kinova/ArmourDiscrete/include)

target_link_libraries(Kinova_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(Kinova_example PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(KinovaWaitr_example 
    Examples/Kinova/WaitrDiscrete/KinovaWaitrExample.cpp)

target_include_directories(KinovaWaitr_example PUBLIC
    Examples/Kinova
    Examples/Kinova/WaitrDiscrete/include)

target_link_libraries(KinovaWaitr_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(KinovaWaitr_example PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(KinovaFilter_example 
    Examples/Kinova/SystemIdentification/DataFilter/KinovaAccelerationEstimator.cpp)

target_include_directories(KinovaFilter_example PUBLIC
    Examples/Kinova
    Examples/Kinova/SystemIdentification/DataFilter/include
    Examples/Kinova/SystemIdentification/ExcitingTrajectories/include
    Examples/Kinova/SystemIdentification/IterativeSystemIdentification/include)

target_link_libraries(KinovaFilter_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES})

add_executable(KinovaSysId_example 
    Examples/Kinova/SystemIdentification/ExcitingTrajectories/KinovaRegressorExample.cpp)

target_include_directories(KinovaSysId_example PUBLIC
    Examples/Kinova
    Examples/Kinova/SystemIdentification/DataFilter/include
    Examples/Kinova/SystemIdentification/ExcitingTrajectories/include
    Examples/Kinova/SystemIdentification/IterativeSystemIdentification/include)

target_link_libraries(KinovaSysId_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(KinovaSysId_example PUBLIC
    ${PINOCCHIO_FLAGS})

nanobind_add_module(oracle_nanobind 
    NB_SHARED LTO
    Examples/Kinova/ArmourDiscrete/KinovaPybind.cpp
    Examples/Kinova/ArmourDiscrete/src/KinovaPybindWrapper.cpp)

target_include_directories(oracle_nanobind PUBLIC
    Examples/Kinova
    Examples/Kinova/ArmourDiscrete/include)

target_link_libraries(oracle_nanobind PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio
    ${PYTHON_LIBRARIES})

nanobind_compile_options(oracle_nanobind PUBLIC
    ${PINOCCHIO_FLAGS})

set_property(TARGET oracle_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)

nanobind_add_module(oracle_waitr_nanobind SHARED
    Examples/Kinova/WaitrDiscrete/KinovaWaitrPybind.cpp
    Examples/Kinova/WaitrDiscrete/src/KinovaWaitrPybindWrapper.cpp)

target_include_directories(oracle_waitr_nanobind PUBLIC
    Examples/Kinova
    Examples/Kinova/WaitrDiscrete/include)

target_link_libraries(oracle_waitr_nanobind PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio
    ${PYTHON_LIBRARIES})

    nanobind_compile_options(oracle_waitr_nanobind PUBLIC
    ${PINOCCHIO_FLAGS})

set_property(TARGET oracle_waitr_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)

## Test scripts
add_executable(ForwardKinematics_test 
    Tests/KinematicsDynamics/TestForwardKinematicsSolver.cpp)

target_link_libraries(ForwardKinematics_test PUBLIC
    IDlib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(ForwardKinematics_test PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(ForwardKinematicsGradient_test 
    Tests/KinematicsDynamics/TestForwardKinematicsRPYDerivatives.cpp)

target_link_libraries(ForwardKinematicsGradient_test PUBLIC
    IDlib
    Optlib
    Conslib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(ForwardKinematicsGradient_test PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(KinematicsConstraints_test 
    Tests/Constraints/TestKinematicsConstraints.cpp)

target_link_libraries(KinematicsConstraints_test PUBLIC
    trajlib
    IDlib
    Conslib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(KinematicsConstraints_test PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(CustomizedInverseDynamics_test 
    Tests/KinematicsDynamics/TestCustomizedInverseDynamics.cpp)

target_link_libraries(CustomizedInverseDynamics_test PUBLIC
    trajlib
    IDlib
    Conslib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(CustomizedInverseDynamics_test PUBLIC
    ${PINOCCHIO_FLAGS})

## Install
install(TARGETS 
    trajlib
    IDlib
    Conslib
    Optlib
    Kinovalib
    Digitlib
    Taloslib
    # DigitModifiedlib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)