set(COMMON_LIBRARIES
    trajlib
    IDlib
    Costslib
    Conslib
    Optlib
    Armourlib
    ${Boost_LIBRARIES}
    pinocchio::pinocchio
    ipopt
    coinhsl
)

add_executable(ForwardKinematics_test 
    KinematicsDynamics/TestForwardKinematicsSolver.cpp)

target_link_libraries(ForwardKinematics_test PUBLIC
    IDlib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(ForwardKinematics_test PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(ForwardKinematicsGradient_test 
    KinematicsDynamics/TestForwardKinematicsRPYDerivatives.cpp)

target_link_libraries(ForwardKinematicsGradient_test PUBLIC
    IDlib
    Optlib
    Conslib
    ipopt
    coinhsl
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(ForwardKinematicsGradient_test PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(MinimizeTorque_test 
    Costs/TestMinimizeTorqueCosts.cpp)

target_link_libraries(MinimizeTorque_test PUBLIC
    trajlib
    IDlib
    Costslib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(MinimizeTorque_test PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(KinematicsConstraints_test 
    Constraints/TestKinematicsConstraints.cpp)

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
    KinematicsDynamics/TestCustomizedInverseDynamics.cpp)

target_link_libraries(CustomizedInverseDynamics_test PUBLIC
    trajlib
    IDlib
    Conslib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(CustomizedInverseDynamics_test PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(RegressorInverseDynamics_test 
    KinematicsDynamics/TestRegressorInverseDynamics.cpp)

target_link_libraries(RegressorInverseDynamics_test PUBLIC
    trajlib
    IDlib
    Conslib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(RegressorInverseDynamics_test PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(PZsparse_test 
    Armour/PZsparseTest.cpp)

target_link_libraries(PZsparse_test PUBLIC
    trajlib
    IDlib
    Conslib
    Optlib
    Armourlib
    ${BOOST_LIBRARIES}
    yaml-cpp
    pinocchio::pinocchio)

add_executable(ReachableSets_test 
    Armour/ReachableSetsTest.cpp)

target_link_libraries(ReachableSets_test PUBLIC
    trajlib
    IDlib
    Conslib
    Optlib
    Armourlib
    ${BOOST_LIBRARIES}
    yaml-cpp
    pinocchio::pinocchio)

add_executable(QRDecomposition_test
    KinematicsDynamics/TestQRDecomposition.cpp)

target_link_libraries(QRDecomposition_test PUBLIC
    trajlib 
    IDlib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(QRDecomposition_test PUBLIC
    ${PINOCCHIO_FLAGS})

# add_executable(Test_torque
#     KinematicsDynamics/Test_torque.cpp)

# target_link_libraries(Test_torque PUBLIC
#     trajlib 
#     IDlib
#     ${BOOST_LIBRARIES} 
#     pinocchio::pinocchio)

# target_compile_options(Test_torque PUBLIC
#     ${PINOCCHIO_FLAGS})