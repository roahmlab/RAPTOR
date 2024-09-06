# Kinova 库
add_library(Kinovalib SHARED
    InverseKinematics/KinovaIKSolver.cpp
    ArmourDiscrete/src/KinovaOptimizer.cpp
    ArmourDiscrete/src/KinovaCustomizedConstraints.cpp
    WaitrDiscrete/src/KinovaWaitrOptimizer.cpp
    SystemIdentification/DataFilter/src/DataFilterOptimizer.cpp
    SystemIdentification/ExcitingTrajectories/src/ConditionNumberOptimizer.cpp
)

target_include_directories(Kinovalib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/InverseKinematics
    ${CMAKE_CURRENT_SOURCE_DIR}/ArmourDiscrete/include
    ${CMAKE_CURRENT_SOURCE_DIR}/WaitrDiscrete/include
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/DataFilter/include
    ${CMAKE_CURRENT_SOURCE_DIR}/SystemIdentification/ExcitingTrajectories/include
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
)

# Kinova IK 
add_executable(KinovaIK_example 
    InverseKinematics/KinovaIKExample.cpp
)

target_link_libraries(KinovaIK_example PUBLIC
    Kinovalib
)

# Kinova Armour 
add_executable(KinovaArmour_example 
    ArmourDiscrete/KinovaExample.cpp
)

target_link_libraries(KinovaArmour_example PUBLIC
    Kinovalib
)

# Kinova Waitr 
add_executable(KinovaWaitr_example 
    WaitrDiscrete/KinovaWaitrExample.cpp
)

target_link_libraries(KinovaWaitr_example PUBLIC
    Kinovalib
)

# Kinova System Identification 
add_executable(KinovaSysId_example 
    SystemIdentification/ExcitingTrajectories/KinovaRegressorExample.cpp
)

target_link_libraries(KinovaSysId_example PUBLIC
    Kinovalib
)

# Python bindings 
if(DEFINED PYTHON_BINDING)
    nanobind_add_module(KinovaHLP_nanobind 
        NB_SHARED LTO
        src/WaypointPlanning/WaypointPlanningPybind.cpp
    )

    target_link_libraries(KinovaHLP_nanobind PUBLIC
        Kinovalib
        ${PYTHON_LIBRARIES}
    )

    set_property(TARGET KinovaHLP_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()