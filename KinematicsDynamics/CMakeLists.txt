add_library(IDlib SHARED
    src/InverseDynamics.cpp
    src/CustomizedInverseDynamics.cpp
    src/RegressorInverseDynamics.cpp
    src/MomentumRegressor.cpp
    src/IntervalRegressorInverseDynamics.cpp
    src/IntervalMomentumRegressor.cpp
    src/DynamicsConstraints.cpp
    src/ConstrainedInverseDynamics.cpp
    src/Transform.cpp
    src/ForwardKinematics.cpp
)

target_include_directories(IDlib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/Trajectories/include
    ${CMAKE_SOURCE_DIR}/Utils/include
)