# Minimum required CMake version
cmake_minimum_required(VERSION 3.18)

# Project name
project(IDTO
    DESCRIPTION "Inverse Dynamics Trajectory Optimization")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp")

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find Python library
add_subdirectory(extern/pybind)

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
# have to manually set the include and library directories for yaml-cpp for now
set(YAML_CPP_INCLUDE_DIR /usr/include)
set(YAML_CPP_LIBRARY_DIR /usr/lib/x86_64-linux-gnu/)
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
    Trajectories/src/WaitrBezierCurves.cpp
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
    Examples/Digit/src/DigitMultipleStepDynamicsConstraints.cpp
    Examples/Digit/src/DigitConstrainedInverseDynamics.cpp
    Examples/Digit/src/DigitSingleStepOptimizer.cpp
    Examples/Digit/src/DigitMultipleStepOptimizer.cpp
    Examples/Digit/src/DigitCustomizedConstraints.cpp
    Examples/Digit/src/DigitMultipleStepCustomizedConstraints.cpp
    Examples/Digit/src/DigitSingleStepPeriodicityConstraints.cpp)

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

## Digit-modified library
add_library(DigitModifiedlib SHARED
    Examples/Digit-modified/src/DigitModifiedDynamicsConstraints.cpp
    Examples/Digit-modified/src/DigitModifiedConstrainedInverseDynamics.cpp
    Examples/Digit-modified/src/DigitModifiedSingleStepOptimizer.cpp
    Examples/Digit-modified/src/DigitModifiedCustomizedConstraints.cpp
    Examples/Digit-modified/src/DigitModifiedSingleStepPeriodicityConstraints.cpp)

target_link_libraries(DigitModifiedlib PUBLIC
    trajlib
    IDlib
    Conslib
    Optlib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_include_directories(DigitModifiedlib PUBLIC
    Examples/Digit-modified/include)

# Kinova libraries
add_library(Kinovalib SHARED
    Examples/Kinova/InverseKinematics/KinovaIKSolver.cpp
    Examples/Kinova/ArmourUnsafe/src/KinovaOptimizer.cpp
    Examples/Kinova/ArmourUnsafe/src/KinovaCustomizedConstraints.cpp
    Examples/Kinova/WaitrUnsafe/src/KinovaWaitrOptimizer.cpp
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
    Examples/Kinova/ArmourUnsafe/include
    Examples/Kinova/WaitrUnsafe/include
    Examples/Kinova/SystemIdentification/DataFilter/include
    Examples/Kinova/SystemIdentification/ExcitingTrajectories/include
    Examples/Kinova/SystemIdentification/IterativeSystemIdentification/include)

## Example scripts
add_executable(Digit_example 
    Examples/Digit/DigitSingleStep.cpp)

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
    ${GSL_LIBRARIES}
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(Digit_example PUBLIC
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
    ${GSL_LIBRARIES}
    ${BOOST_LIBRARIES} 
    ${PKG_PIN_CONFIG_LIBRARIES})

target_compile_options(DigitMultipleStep_example PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(DigitModified_example 
    Examples/Digit-modified/DigitModifiedSingleStep.cpp)

target_include_directories(DigitModified_example PUBLIC
    Examples/Digit-modified/include)

target_link_libraries(DigitModified_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    DigitModifiedlib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(DigitModified_example PUBLIC
    ${PINOCCHIO_FLAGS})

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
    Examples/Kinova/ArmourUnsafe/KinovaExample.cpp)

target_include_directories(Kinova_example PUBLIC
    Examples/Kinova
    Examples/Kinova/ArmourUnsafe/include)

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
    Examples/Kinova/WaitrUnsafe/KinovaWaitrExample.cpp)

target_include_directories(KinovaWaitr_example PUBLIC
    Examples/Kinova
    Examples/Kinova/WaitrUnsafe/include)

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

pybind11_add_module(oracle_pybind SHARED
    Examples/Kinova/ArmourUnsafe/KinovaPybind.cpp
    Examples/Kinova/ArmourUnsafe/src/KinovaPybindWrapper.cpp)

target_include_directories(oracle_pybind PUBLIC
    Examples/Kinova
    Examples/Kinova/ArmourUnsafe/include)

target_link_libraries(oracle_pybind PUBLIC
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

target_compile_options(oracle_pybind PUBLIC
    ${PINOCCHIO_FLAGS})

set_property(TARGET oracle_pybind PROPERTY POSITION_INDEPENDENT_CODE ON)

pybind11_add_module(oracle_waitr_pybind SHARED
    Examples/Kinova/WaitrUnsafe/KinovaWaitrPybind.cpp
    Examples/Kinova/WaitrUnsafe/src/KinovaWaitrPybindWrapper.cpp)

target_include_directories(oracle_waitr_pybind PUBLIC
    Examples/Kinova
    Examples/Kinova/WaitrUnsafe/include)

target_link_libraries(oracle_waitr_pybind PUBLIC
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

target_compile_options(oracle_waitr_pybind PUBLIC
    ${PINOCCHIO_FLAGS})

set_property(TARGET oracle_waitr_pybind PROPERTY POSITION_INDEPENDENT_CODE ON)

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

## Install
install(TARGETS 
    trajlib
    IDlib
    Conslib
    Optlib
    Kinovalib
    Digitlib
    DigitModifiedlib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)