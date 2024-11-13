add_library(IDlib SHARED
    src/InverseDynamics.cpp
    src/CustomizedInverseDynamics.cpp
    src/QRDecompositionSolver.cpp
    src/RegressorInverseDynamics.cpp
    src/MomentumRegressor.cpp
    src/DynamicsConstraints.cpp
    src/ConstrainedInverseDynamics.cpp
    src/Transform.cpp
    src/Spatial.cpp
    src/ForwardKinematics.cpp
)

target_include_directories(IDlib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/Trajectories/include
    ${CMAKE_SOURCE_DIR}/Utils/include
)