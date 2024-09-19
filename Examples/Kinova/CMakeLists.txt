# Kinova 
add_library(Kinovalib SHARED
    InverseKinematics/KinovaIKSolver.cpp
    ArmourDiscrete/src/KinovaOptimizer.cpp
    ArmourDiscrete/src/KinovaCustomizedConstraints.cpp
    WaitrDiscrete/src/KinovaWaitrOptimizer.cpp
    SystemIdentification/DataFilter/src/DataFilterOptimizer.cpp
    SystemIdentification/ExcitingTrajectories/src/ConditionNumberOptimizer.cpp
    SystemIdentification/ExcitingTrajectories/src/EndEffectorConditionalNumberOptimizer.cpp
)

target_include_directories(Kinovalib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/InverseKinematics
    ${CMAKE_CURRENT_SOURCE_DIR}/ArmourDiscrete/include
    ${CMAKE_CURRENT_SOURCE_DIR}/WaitrDiscrete/include
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/DataFilter/include
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/ExcitingTrajectories/include
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/IterativeSystemIdentification/include
)

target_link_libraries(Kinovalib PUBLIC
    trajlib
    IDlib
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
    Armour/src/PZsparse.cpp
    Armour/src/RobotInfo.cpp
    Armour/src/ParameterizedTrajectory.cpp
    Armour/src/PZDynamics.cpp
    Armour/src/ReachableSets.cpp
)

target_include_directories(Armourlib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Armour/include
)

target_link_libraries(Armourlib PUBLIC
    ${Boost_LIBRARIES} 
    yaml-cpp
    pinocchio::pinocchio
)

###
add_executable(KinovaIK_example 
    InverseKinematics/KinovaIKExample.cpp)

target_include_directories(KinovaIK_example PUBLIC
    InverseKinematics)

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
add_executable(KinovaHLP_example 
    WaypointPlanning/WaypointPlanningExample.cpp)

target_include_directories(KinovaHLP_example PUBLIC
    ArmourDiscrete/include)

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
    ArmourDiscrete/KinovaExample.cpp)

target_include_directories(Kinova_example PUBLIC
    ArmourDiscrete/include)

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
add_executable(KinovaFilter_example 
    SystemIdentification/DataFilter/KinovaAccelerationEstimator.cpp)


target_include_directories(KinovaFilter_example PUBLIC
    systemIdentification/DataFilter/include
    SystemIdentification/ExcitingTrajectories/include
    SystemIdentification/IterativeSystemIdentification/include)

target_link_libraries(KinovaFilter_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES})


###
add_executable(KinovaExciting_example 
    SystemIdentification/ExcitingTrajectories/KinovaRegressorExample.cpp)

target_include_directories(KinovaExciting_example PUBLIC
    ystemIdentification/DataFilter/include
    SystemIdentification/ExcitingTrajectories/include
    SystemIdentification/IterativeSystemIdentification/include)

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


add_executable(KinovaEndeffectorExciting_example
    SystemIdentification/ExcitingTrajectories/KinovaRegressorExampleEndeffector.cpp)

target_include_directories(KinovaEndeffectorExciting_example PUBLIC
    ystemIdentification/DataFilter/include
    SystemIdentification/ExcitingTrajectories/include
    SystemIdentification/IterativeSystemIdentification/include)

target_link_libraries(KinovaEndeffectorExciting_example PUBLIC
    trajlib 
    IDlib
    Conslib
    Optlib
    Kinovalib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(KinovaEndeffectorExciting_example PUBLIC
    ${PINOCCHIO_FLAGS})


### Python bindings 
nanobind_add_module(KinovaHLP_nanobind 
    NB_SHARED LTO
    WaypointPlanning/WaypointPlanningPybind.cpp)

target_include_directories(KinovaHLP_nanobind PUBLIC
    ArmourDiscrete/include)

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

nanobind_add_module(oracle_nanobind 
    NB_SHARED LTO
    ArmourDiscrete/KinovaPybind.cpp
    ArmourDiscrete/src/KinovaPybindWrapper.cpp)

target_include_directories(oracle_nanobind PUBLIC
    Kinova/ArmourDiscrete/include)

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
    WaitrDiscrete/KinovaWaitrPybind.cpp
    WaitrDiscrete/src/KinovaWaitrPybindWrapper.cpp)

target_include_directories(oracle_waitr_nanobind PUBLIC
    WaitrDiscrete/include)

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