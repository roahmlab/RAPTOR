# Kinova 
add_library(Kinovalib SHARED
    CollisionAvoidanceInverseKinematics/src/KinovaIKSolver.cpp
    CollisionAvoidanceTrajectory/src/KinovaOptimizer.cpp
    CollisionAvoidanceTrajectory/src/KinovaCustomizedConstraints.cpp
    CollisionAvoidanceTrajectory/src/KinovaLongerHorizonOptimizer.cpp
    CollisionAvoidanceTrajectory/src/DualKinovaOptimizer.cpp
    WaitrDiscrete/src/KinovaWaitrOptimizer.cpp
    SystemIdentification/ExcitingTrajectories/src/ExcitingTrajectoryGenerator.cpp
    SystemIdentification/ExcitingTrajectories/src/PayloadExcitingTrajectoryGenerator.cpp
    SystemIdentification/ParametersIdentification/src/FrictionParametersIdentification.cpp
    SystemIdentification/ParametersIdentification/src/EndEffectorParametersIdentification.cpp
    SystemIdentification/ParametersIdentification/src/EndEffectorParametersIdentificationMomentum.cpp)
    
target_include_directories(Kinovalib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/CollisionAvoidanceInverseKinematics/include
    ${CMAKE_CURRENT_SOURCE_DIR}/CollisionAvoidanceTrajectory/include
    ${CMAKE_CURRENT_SOURCE_DIR}/WaitrDiscrete/include
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/ExcitingTrajectories/include
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/ParametersIdentification/include
    ${CMAKE_SOURCE_DIR}/Trajectories/include
    ${CMAKE_SOURCE_DIR}/PZDynamics/include
    ${CMAKE_SOURCE_DIR}/Costs/include
    ${CMAKE_SOURCE_DIR}/Constraints/include
    ${CMAKE_SOURCE_DIR}/Optimization/include
    ${CMAKE_SOURCE_DIR}/Utils/include
)

target_link_libraries(Kinovalib PUBLIC
    trajlib
    IDlib
    Costslib
    Conslib
    Optlib
    ipopt
    coinhsl
    ${Boost_LIBRARIES} 
    pinocchio::pinocchio
    ompl
)

# Armour
add_library(Armourlib SHARED
    Armour/src/PZSparse.cpp
    Armour/src/RobotInfo.cpp
    Armour/src/ParameterizedTrajectory.cpp
    Armour/src/PZDynamics.cpp
    Armour/src/ArmourOptimizer.cpp
    Armour/src/DualArmourOptimizer.cpp
    SystemIdentification/ExcitingTrajectories/src/SafePayloadExcitingTrajectoryGenerator.cpp
)

target_include_directories(Armourlib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Armour/include
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/ExcitingTrajectories/include
    ${CMAKE_SOURCE_DIR}/Trajectories/include
    ${CMAKE_SOURCE_DIR}/PZDynamics/include
    ${CMAKE_SOURCE_DIR}/Constraints/include
    ${CMAKE_SOURCE_DIR}/Optimization/include
    ${CMAKE_SOURCE_DIR}/Utils/include
)

target_link_libraries(Armourlib PUBLIC
    trajlib
    IDlib
    Conslib
    Costslib
    Optlib
    ipopt
    coinhsl
    ${Boost_LIBRARIES} 
    pinocchio::pinocchio
    yaml-cpp
)

###
add_executable(KinovaIK_example 
    CollisionAvoidanceInverseKinematics/KinovaIKExample.cpp)
target_include_directories(KinovaIK_example PUBLIC
    CollisionAvoidanceInverseKinematics)
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


###
add_executable(KinovaIKMotion_example 
    CollisionAvoidanceInverseKinematics/KinovaIKMotionExample.cpp)
target_include_directories(KinovaIKMotion_example PUBLIC
    CollisionAvoidanceInverseKinematics/include)
target_link_libraries(KinovaIKMotion_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)
target_compile_options(KinovaIKMotion_example PUBLIC
    ${PINOCCHIO_FLAGS})


###
add_executable(KinovaHLP_example 
    WaypointPlanning/WaypointPlanningExample.cpp)
target_include_directories(KinovaHLP_example PUBLIC
    CollisionAvoidanceTrajectory/include)
target_link_libraries(KinovaHLP_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Kinovalib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio
    ompl)
target_compile_options(KinovaHLP_example PUBLIC
    ${PINOCCHIO_FLAGS})


###
add_executable(Kinova_example 
    CollisionAvoidanceTrajectory/KinovaExample.cpp)
target_include_directories(Kinova_example PUBLIC
    CollisionAvoidanceTrajectory/include)
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


###
add_executable(KinovaWaitr_example 
    WaitrDiscrete/KinovaWaitrExample.cpp)
target_include_directories(KinovaWaitr_example PUBLIC
    WaitrDiscrete/include)
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


###
add_executable(KinovaExciting_example 
    SystemIdentification/ExcitingTrajectories/KinovaRegressorExample.cpp)
target_include_directories(KinovaExciting_example PUBLIC
    SystemIdentification/ExcitingTrajectories/include
    SystemIdentification/ParametersIdentification/include)
target_link_libraries(KinovaExciting_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)
target_compile_options(KinovaExciting_example PUBLIC
    ${PINOCCHIO_FLAGS})


###
add_executable(KinovaPayloadExciting_example
    SystemIdentification/ExcitingTrajectories/KinovaRegressorExamplePayload.cpp)
target_include_directories(KinovaPayloadExciting_example PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/ExcitingTrajectories/include
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/ParametersIdentification/include)
target_link_libraries(KinovaPayloadExciting_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)
target_compile_options(KinovaPayloadExciting_example PUBLIC
    ${PINOCCHIO_FLAGS})


###
add_executable(DualKinova_example 
    CollisionAvoidanceTrajectory/DualKinovaExample.cpp)
target_include_directories(DualKinova_example PUBLIC
    CollisionAvoidanceTrajectory/include)
target_link_libraries(DualKinova_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Costslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)
target_compile_options(DualKinova_example PUBLIC
    ${PINOCCHIO_FLAGS})


###
add_executable(SysidFriction_test
    SystemIdentification/ParametersIdentification/TestFrictionParametersIdentification.cpp)
target_include_directories(SysidFriction_test PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/ParametersIdentification/include)
target_link_libraries(SysidFriction_test PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)
target_compile_options(SysidFriction_test PUBLIC
    ${PINOCCHIO_FLAGS})


###
add_executable(SysidEnd_test
    SystemIdentification/ParametersIdentification/TestEndEffectorParametersIdentification.cpp)
target_include_directories(SysidEnd_test PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/ParametersIdentification/include)
target_link_libraries(SysidEnd_test PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)
target_compile_options(SysidEnd_test PUBLIC
    ${PINOCCHIO_FLAGS})


###
add_executable(SysidEndMomentum_test
    SystemIdentification/ParametersIdentification/TestEndEffectorParametersIdentificationMomentum.cpp)
target_include_directories(SysidEndMomentum_test PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/ParametersIdentification/include)
target_link_libraries(SysidEndMomentum_test PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)
target_compile_options(SysidEndMomentum_test PUBLIC
    ${PINOCCHIO_FLAGS})
    

###
add_executable(Armour_example 
    Armour/ArmourExample.cpp)
target_include_directories(Armour_example PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/Armour/include)
target_link_libraries(Armour_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Armourlib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)
target_compile_options(Armour_example PUBLIC
    ${PINOCCHIO_FLAGS})


add_executable(DualArmour_example 
    Armour/DualArmourExample.cpp)
target_include_directories(DualArmour_example PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/Armour/include)
target_link_libraries(DualArmour_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Armourlib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)
target_compile_options(DualArmour_example PUBLIC
    ${PINOCCHIO_FLAGS})


### Python bindings 
nanobind_add_module(KinovaIKMotion_nanobind 
    NB_SHARED LTO
    CollisionAvoidanceInverseKinematics/KinovaIKMotionPybind.cpp
    CollisionAvoidanceInverseKinematics/src/KinovaIKMotionPybindWrapper.cpp)
target_include_directories(KinovaIKMotion_nanobind PUBLIC
    CollisionAvoidanceInverseKinematics/include)
target_link_libraries(KinovaIKMotion_nanobind PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio
    ${PYTHON_LIBRARIES})
nanobind_compile_options(KinovaIKMotion_nanobind PUBLIC
    ${PINOCCHIO_FLAGS})
set_property(TARGET KinovaIKMotion_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)


nanobind_add_module(KinovaHLP_nanobind 
    NB_SHARED LTO
    WaypointPlanning/WaypointPlanningPybind.cpp)
target_include_directories(KinovaHLP_nanobind PUBLIC
    CollisionAvoidanceTrajectory/include)
target_link_libraries(KinovaHLP_nanobind PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio
    ompl
    ${PYTHON_LIBRARIES})
nanobind_compile_options(KinovaHLP_nanobind PUBLIC
    ${PINOCCHIO_FLAGS})
set_property(TARGET KinovaHLP_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)


nanobind_add_module(armour_nanobind 
    NB_SHARED LTO
    Armour/ArmourPybind.cpp
    Armour/src/ArmourPybindWrapper.cpp)
target_include_directories(armour_nanobind PUBLIC
    Armour/include)
target_link_libraries(armour_nanobind PUBLIC
    # trajlib 
    # IDlib
    # Conslib
    # Optlib
    Kinovalib
    Armourlib
    # ipopt
    # coinhsl
    # ${BOOST_LIBRARIES} 
    # pinocchio::pinocchio
    ${PYTHON_LIBRARIES})
nanobind_compile_options(armour_nanobind PUBLIC
    ${PINOCCHIO_FLAGS})
set_property(TARGET armour_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)


nanobind_add_module(dualarmour_nanobind 
    NB_SHARED LTO
    Armour/DualArmourPybind.cpp
    Armour/src/DualArmourPybindWrapper.cpp)
target_include_directories(dualarmour_nanobind PUBLIC
    Armour/include)
target_link_libraries(dualarmour_nanobind PUBLIC
    # trajlib 
    # IDlib
    # Conslib
    # Optlib
    Kinovalib
    Armourlib
    # ipopt
    # coinhsl
    # ${BOOST_LIBRARIES} 
    # pinocchio::pinocchio
    ${PYTHON_LIBRARIES})
nanobind_compile_options(dualarmour_nanobind PUBLIC
    ${PINOCCHIO_FLAGS})
set_property(TARGET dualarmour_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)


nanobind_add_module(kinova_nanobind 
    NB_SHARED LTO
    CollisionAvoidanceTrajectory/KinovaPybind.cpp
    CollisionAvoidanceTrajectory/src/KinovaPybindWrapper.cpp)
target_include_directories(kinova_nanobind PUBLIC
    CollisionAvoidanceTrajectory/include)
target_link_libraries(kinova_nanobind PUBLIC
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
nanobind_compile_options(kinova_nanobind PUBLIC
    ${PINOCCHIO_FLAGS})
set_property(TARGET kinova_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)


nanobind_add_module(kinova_longer_nanobind 
    NB_SHARED LTO
    CollisionAvoidanceTrajectory/KinovaLongerHorizonPybind.cpp
    CollisionAvoidanceTrajectory/src/KinovaLongerHorizonPybindWrapper.cpp)
target_include_directories(kinova_longer_nanobind PUBLIC
    CollisionAvoidanceTrajectory/include)
target_link_libraries(kinova_longer_nanobind PUBLIC
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
nanobind_compile_options(kinova_longer_nanobind PUBLIC
    ${PINOCCHIO_FLAGS})
set_property(TARGET kinova_longer_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)


nanobind_add_module(kinova_waitr_nanobind SHARED
    WaitrDiscrete/KinovaWaitrPybind.cpp
    WaitrDiscrete/src/KinovaWaitrPybindWrapper.cpp)
target_include_directories(kinova_waitr_nanobind PUBLIC
    WaitrDiscrete/include)
target_link_libraries(kinova_waitr_nanobind PUBLIC
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
nanobind_compile_options(kinova_waitr_nanobind PUBLIC
    ${PINOCCHIO_FLAGS})
set_property(TARGET kinova_waitr_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)


nanobind_add_module(safe_payload_exciting_nanobind 
    NB_SHARED LTO
    SystemIdentification/ExcitingTrajectories/SafePayloadExcitingPybind.cpp
    SystemIdentification/ExcitingTrajectories/src/SafePayloadExcitingPybindWrapper.cpp
    SystemIdentification/ExcitingTrajectories/src/SafePayloadExcitingTrajectoryGenerator.cpp)
target_include_directories(safe_payload_exciting_nanobind PUBLIC
    Armour/include
    SystemIdentification/ExcitingTrajectories/include)
target_link_libraries(safe_payload_exciting_nanobind PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    Armourlib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio
    ${PYTHON_LIBRARIES})
nanobind_compile_options(safe_payload_exciting_nanobind PUBLIC
    ${PINOCCHIO_FLAGS})
set_property(TARGET safe_payload_exciting_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)


nanobind_add_module(end_effector_sysid_nanobind 
    NB_SHARED LTO
    SystemIdentification/ParametersIdentification/src/EndEffectorIdentificationPybindWrapper.cpp
    SystemIdentification/ParametersIdentification/EndEffectorIdentificationPybind.cpp)
target_include_directories(end_effector_sysid_nanobind PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/ParametersIdentification/include)
target_link_libraries(end_effector_sysid_nanobind PUBLIC
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
nanobind_compile_options(end_effector_sysid_nanobind PUBLIC
    ${PINOCCHIO_FLAGS})
set_property(TARGET end_effector_sysid_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)


nanobind_add_module(end_effector_sysid_momentum_nanobind 
    NB_SHARED LTO
    SystemIdentification/ParametersIdentification/src/EndEffectorIdentificationMomentumPybindWrapper.cpp
    SystemIdentification/ParametersIdentification/EndEffectorIdentificationMomentumPybind.cpp)
target_include_directories(end_effector_sysid_momentum_nanobind PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/ParametersIdentification/include)
target_link_libraries(end_effector_sysid_momentum_nanobind PUBLIC
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
nanobind_compile_options(end_effector_sysid_momentum_nanobind PUBLIC
    ${PINOCCHIO_FLAGS})
set_property(TARGET end_effector_sysid_momentum_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)
