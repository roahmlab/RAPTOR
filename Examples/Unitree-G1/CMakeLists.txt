add_library(G1lib SHARED
    src/G1DynamicsConstraints.cpp
    src/G1ConstrainedInverseDynamics.cpp
    src/G1CustomizedConstraints.cpp
    src/G1SingleStepOptimizer.cpp
    src/G1SingleStepPeriodicityConstraints.cpp
    src/G1MultipleStepOptimizer.cpp
    src/G1MultipleStepPeriodicityConstraints.cpp
)

target_include_directories(G1lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/KinematicsDynamics/include
    ${CMAKE_SOURCE_DIR}/Trajectories/include
    ${CMAKE_SOURCE_DIR}/Constraints/include
    ${CMAKE_SOURCE_DIR}/Costs/include
    ${CMAKE_SOURCE_DIR}/Optimization/include
    ${CMAKE_SOURCE_DIR}/Utils/include)

target_link_libraries(G1lib PUBLIC
    trajlib
    IDlib
    Conslib
    Costslib
    Optlib
    ipopt
    coinhsl
    ${Boost_LIBRARIES} 
    pinocchio::pinocchio
)

add_executable(G1_example 
    G1SingleStep.cpp)
    # VisualizationHelper.cpp)

target_link_libraries(G1_example PUBLIC
    G1lib
    yaml-cpp
)
target_compile_options(G1_example PUBLIC
    ${PINOCCHIO_FLAGS})

# add_executable(G1MultipleStep_example 
#     G1MultipleStep.cpp
# )

# target_link_libraries(G1MultipleStep_example PUBLIC
#     G1lib
#     yaml-cpp
# )
# target_compile_options(G1MultipleStep_example PUBLIC
#     ${PINOCCHIO_FLAGS})