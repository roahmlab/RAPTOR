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
)

target_link_libraries(Taloslib PUBLIC
    trajlib
    IDlib
    Conslib
    Optlib
    ipopt
    coinhsl
    ${Boost_LIBRARIES} 
    pinocchio::pinocchio
)

#execute
add_executable(TalosSingleStep_example 
    TalosSingleStep.cpp
)

target_link_libraries(TalosSingleStep_example PUBLIC
    Taloslib
    yaml-cpp
)


add_executable(TalosMultipleStep_example 
    TalosMultipleStep.cpp
)

target_link_libraries(TalosMultipleStep_example PUBLIC
    Taloslib
    yaml-cpp
)


add_executable(TalosSingleStepFixedPosition_example 
    TalosSingleStepFixedPosition.cpp
)

target_link_libraries(TalosSingleStepFixedPosition_example PUBLIC
    Taloslib
    yaml-cpp
)