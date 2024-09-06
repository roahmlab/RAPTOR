add_library(Digitlib SHARED
    src/DigitDynamicsConstraints.cpp
    src/DigitConstrainedInverseDynamics.cpp
    src/DigitSingleStepOptimizer.cpp
    src/DigitMultipleStepOptimizer.cpp
    src/DigitCustomizedConstraints.cpp
    src/DigitSingleStepPeriodicityConstraints.cpp
    src/DigitMultipleStepPeriodicityConstraints.cpp
)

target_include_directories(Digitlib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(Digitlib PUBLIC
    trajlib
    IDlib
    Conslib
    Optlib
    ipopt
    coinhsl
    ${GSL_LIBRARIES}
    ${Boost_LIBRARIES} 
    pinocchio::pinocchio
)

add_executable(Digit_example 
    DigitSingleStep.cpp
)

target_link_libraries(Digit_example PUBLIC
    Digitlib
    yaml-cpp
)