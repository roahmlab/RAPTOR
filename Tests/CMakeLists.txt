# Define common libraries
set(COMMON_LIBRARIES
    trajlib
    IDlib
    Conslib
    Optlib
    ${Boost_LIBRARIES}
    pinocchio::pinocchio
    ipopt
    coinhsl
)

# Define a function to add tests
function(add_raptor_test test_name source_file)
    add_executable(${test_name} ${source_file})
    target_link_libraries(${test_name} PRIVATE ${COMMON_LIBRARIES})
    target_compile_options(${test_name} PRIVATE ${PINOCCHIO_FLAGS})
    add_test(NAME ${test_name} COMMAND ${test_name})
    set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
endfunction()

# Add individual tests
add_raptor_test(ForwardKinematics_test 
                KinematicsDynamics/TestForwardKinematicsSolver.cpp)

add_raptor_test(ForwardKinematicsGradient_test 
        KinematicsDynamics/TestForwardKinematicsRPYDerivatives.cpp)

add_raptor_test(KinematicsConstraints_test 
        Constraints/TestKinematicsConstraints.cpp)

add_raptor_test(CustomizedInverseDynamics_test 
        KinematicsDynamics/TestCustomizedInverseDynamics.cpp)
        
add_raptor_test(RegressorInverseDynamics_test 
            KinematicsDynamics/TestRegressorInverseDynamics.cpp)
