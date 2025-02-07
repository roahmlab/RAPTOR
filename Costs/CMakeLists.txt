add_library(Costslib SHARED
    src/Costs.cpp
    src/MinimizeInitialVelocity.cpp
    src/MinimizeInitialAcceleration.cpp
    src/MinimizeTorque.cpp
    src/MinimizePower.cpp
    src/MinimizePathLength.cpp
    src/MinimizeJerk.cpp
    src/EndEffectorRegressorConditionNumber.cpp
)

target_include_directories(Costslib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/Trajectories/include
    ${CMAKE_SOURCE_DIR}/KinematicsDynamics/include
    ${CMAKE_SOURCE_DIR}/Utils/include
)

target_link_libraries(Costslib PUBLIC
    trajlib
    IDlib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio
)