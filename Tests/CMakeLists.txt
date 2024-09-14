

# set(COMMON_LIBRARIES
#     trajlib
#     IDlib
#     Conslib
#     Optlib
#     ${Boost_LIBRARIES}
#     pinocchio::pinocchio
#     ipopt
#     coinhsl
# )

# add_executable(ForwardKinematics_test 
#     KinematicsDynamics/TestForwardKinematicsSolver.cpp)

# target_link_libraries(ForwardKinematics_test PUBLIC
#     IDlib
#     ${BOOST_LIBRARIES} 
#     pinocchio::pinocchio)

# target_compile_options(ForwardKinematics_test PUBLIC
#     ${PINOCCHIO_FLAGS})

# add_test(NAME ForwardKinematics_test COMMAND ForwardKinematics_test)
# set_tests_properties(ForwardKinematics_test PROPERTIES WORKING_DIRECTORY /workspaces/RAPTOR/build)

# add_executable(ForwardKinematicsGradient_test 
#     KinematicsDynamics/TestForwardKinematicsRPYDerivatives.cpp)

# target_link_libraries(ForwardKinematicsGradient_test PUBLIC
#     IDlib
#     Optlib
#     Conslib
#     ipopt
#     coinhsl
#     ${BOOST_LIBRARIES} 
#     pinocchio::pinocchio)

# target_compile_options(ForwardKinematicsGradient_test PUBLIC
#     ${PINOCCHIO_FLAGS})

# add_test(NAME ForwardKinematicsGradient_test COMMAND ForwardKinematicsGradient_test)
# set_tests_properties(ForwardKinematicsGradient_test PROPERTIES WORKING_DIRECTORY /workspaces/RAPTOR/build)


# add_executable(KinematicsConstraints_test 
#     Constraints/TestKinematicsConstraints.cpp)

# target_link_libraries(KinematicsConstraints_test PUBLIC
#     trajlib
#     IDlib
#     Conslib
#     ${BOOST_LIBRARIES} 
#     pinocchio::pinocchio)

# target_compile_options(KinematicsConstraints_test PUBLIC
#     ${PINOCCHIO_FLAGS})

# add_test(NAME KinematicsConstraints_test COMMAND KinematicsConstraints_test)
# set_tests_properties(KinematicsConstraints_test PROPERTIES WORKING_DIRECTORY /workspaces/RAPTOR/build)

# add_executable(CustomizedInverseDynamics_test 
#     KinematicsDynamics/TestCustomizedInverseDynamics.cpp)

# target_link_libraries(CustomizedInverseDynamics_test PUBLIC
#     trajlib
#     IDlib
#     Conslib
#     ${BOOST_LIBRARIES} 
#     pinocchio::pinocchio)

# target_compile_options(CustomizedInverseDynamics_test PUBLIC
#     ${PINOCCHIO_FLAGS})

# add_test(NAME CustomizedInverseDynamics_test COMMAND CustomizedInverseDynamics_test)
# set_tests_properties(CustomizedInverseDynamics_test PROPERTIES WORKING_DIRECTORY /workspaces/RAPTOR/build)

# add_executable(RegressorInverseDynamics_test 
#     KinematicsDynamics/TestRegressorInverseDynamics.cpp)

# target_link_libraries(RegressorInverseDynamics_test PUBLIC
#     trajlib
#     IDlib
#     Conslib
#     ${BOOST_LIBRARIES} 
#     pinocchio::pinocchio)

# target_compile_options(RegressorInverseDynamics_test PUBLIC
#     ${PINOCCHIO_FLAGS})
# add_test(NAME RegressorInverseDynamics_test COMMAND RegressorInverseDynamics_test)
# set_tests_properties(RegressorInverseDynamics_test PROPERTIES WORKING_DIRECTORY /workspaces/RAPTOR/build)



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
