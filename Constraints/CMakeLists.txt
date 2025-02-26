add_library(Conslib SHARED
    src/Constraints.cpp
    src/TrajectoryTerminalConstraints.cpp
    src/JointLimits.cpp
    src/VelocityLimits.cpp
    src/AccelerationLimits.cpp
    src/ConstrainedJointLimits.cpp
    src/TorqueLimits.cpp
    src/RectangleSurfaceContactConstraints.cpp
    src/CircleSurfaceContactConstraints.cpp
    src/LieSpaceResidual.cpp
    src/KinematicsConstraints.cpp
    src/BoxCollisionAvoidance.cpp
    src/TaperedCapsuleCollision.cpp
)

target_include_directories(Conslib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/Trajectories/include
    ${CMAKE_SOURCE_DIR}/KinematicsDynamics/include
    ${CMAKE_SOURCE_DIR}/Utils/include
)

target_link_libraries(Conslib PUBLIC
    trajlib
    IDlib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio
)