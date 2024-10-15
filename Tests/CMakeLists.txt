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

<<<<<<< HEAD
add_raptor_test(CustomizedInverseDynamics_test 
=======
target_link_libraries(KinematicsConstraints_test PUBLIC
    trajlib
    IDlib
    Conslib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(KinematicsConstraints_test PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(LMIConstraints_test 
    Constraints/TestLMIConstraints.cpp)

target_link_libraries(LMIConstraints_test PUBLIC
    trajlib
    IDlib
    Conslib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(LMIConstraints_test PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(RegroupedLMIConstraints_test 
    Constraints/TestRegroupedLMIConstraints.cpp)

target_link_libraries(RegroupedLMIConstraints_test PUBLIC
    trajlib
    IDlib
    Conslib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(RegroupedLMIConstraints_test PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(CustomizedInverseDynamics_test 
>>>>>>> 08c148ee62246521b1469efeef65d0c9f39b3b7c
    KinematicsDynamics/TestCustomizedInverseDynamics.cpp)
        
add_raptor_test(RegressorInverseDynamics_test 
    KinematicsDynamics/TestRegressorInverseDynamics.cpp)
<<<<<<< HEAD
=======

target_link_libraries(RegressorInverseDynamics_test PUBLIC
    trajlib
    IDlib
    Conslib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(RegressorInverseDynamics_test PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(QRDecomposition_test
    KinematicsDynamics/TestQRDecomposition.cpp)

target_link_libraries(QRDecomposition_test PUBLIC
    trajlib 
    IDlib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(QRDecomposition_test PUBLIC
    ${PINOCCHIO_FLAGS})
>>>>>>> 08c148ee62246521b1469efeef65d0c9f39b3b7c
