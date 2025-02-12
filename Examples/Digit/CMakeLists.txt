add_library(Digitlib SHARED
    src/DigitDynamicsConstraints.cpp
    src/DigitConstrainedInverseDynamics.cpp
    src/DigitSingleStepOptimizer.cpp
    src/DigitMultipleStepOptimizer.cpp
    src/DigitCustomizedConstraints.cpp
    src/DigitSingleStepPeriodicityConstraints.cpp
    src/DigitMultipleStepPeriodicityConstraints.cpp
    src/DigitWholeBodyDynamicsConstraints.cpp
)

target_include_directories(Digitlib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/KinematicsDynamics/include
    ${CMAKE_SOURCE_DIR}/Trajectories/include
    ${CMAKE_SOURCE_DIR}/Constraints/include
    ${CMAKE_SOURCE_DIR}/Costs/include
    ${CMAKE_SOURCE_DIR}/Optimization/include
    ${CMAKE_SOURCE_DIR}/Utils/include)

target_link_libraries(Digitlib PUBLIC
    trajlib
    IDlib
    Conslib
    Costslib
    Optlib
    ipopt
    coinhsl
    ${GSL_LIBRARIES}
    ${Boost_LIBRARIES} 
    pinocchio::pinocchio
)

add_executable(Digit_example 
    DigitSingleStep.cpp)
    # Examples/Digit/DigitSingleStepSpeedTest.cpp)
    # Examples/Digit/DigitSingleStepRobustnessTest.cpp)
    # Examples/Digit/DigitSingleStepFixedPosition.cpp)
    # Examples/Digit/VisualizationHelper.cpp)

target_link_libraries(Digit_example PUBLIC
    Digitlib
    yaml-cpp
)
target_compile_options(Digit_example PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(DigitMultipleStep_example 
    DigitMultipleStep.cpp)

target_link_libraries(DigitMultipleStep_example PRIVATE
    Digitlib
    yaml-cpp)

target_compile_options(DigitMultipleStep_example PUBLIC
    ${PINOCCHIO_FLAGS})


### Python bindings 
nanobind_add_module(DigitFillUnactuated_nanobind 
    NB_SHARED LTO
    DigitFillUnactuatedJointsPybind.cpp)
target_include_directories(DigitFillUnactuated_nanobind PUBLIC
    include)
target_link_libraries(DigitFillUnactuated_nanobind PUBLIC
    trajlib 
    IDlib
    Conslib
    Costslib
    Digitlib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio
    ${PYTHON_LIBRARIES})
nanobind_compile_options(DigitFillUnactuated_nanobind PUBLIC
    ${PINOCCHIO_FLAGS})
set_property(TARGET DigitFillUnactuated_nanobind PROPERTY POSITION_INDEPENDENT_CODE ON)