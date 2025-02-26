add_library(Taloslib SHARED
    src/TalosDynamicsConstraints.cpp
    src/TalosConstrainedInverseDynamics.cpp
    src/TalosCustomizedConstraints.cpp
    src/TalosSingleStepOptimizer.cpp
    src/TalosSingleStepPeriodicityConstraints.cpp
    src/TalosMultipleStepOptimizer.cpp
    src/TalosMultipleStepPeriodicityConstraints.cpp
)

target_include_directories(Taloslib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/KinematicsDynamics/include
    ${CMAKE_SOURCE_DIR}/Trajectories/include
    ${CMAKE_SOURCE_DIR}/Constraints/include
    ${CMAKE_SOURCE_DIR}/Costs/include
    ${CMAKE_SOURCE_DIR}/Optimization/include
    ${CMAKE_SOURCE_DIR}/Utils/include)

target_link_libraries(Taloslib PUBLIC
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

add_executable(Talos_example 
    TalosSingleStep.cpp)
    # TalosSingleStepFixedPosition.cpp)
    # VisualizationHelper.cpp)

target_link_libraries(Talos_example PUBLIC
    Taloslib
    yaml-cpp
)
target_compile_options(Talos_example PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(TalosMultipleStep_example 
    TalosMultipleStep.cpp
)

target_link_libraries(TalosMultipleStep_example PUBLIC
    Taloslib
    yaml-cpp
)
target_compile_options(TalosMultipleStep_example PUBLIC
    ${PINOCCHIO_FLAGS})