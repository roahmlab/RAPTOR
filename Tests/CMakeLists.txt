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

add_executable(BimanualTest_test 
    Constraints/BimanualTest.cpp)

target_link_libraries(BimanualTest_test PUBLIC
    Armourlib
    ${BOOST_LIBRARIES})

add_executable(BimanualTest_speed 
    Constraints/BimanualTestSpeed.cpp)

target_link_libraries(BimanualTest_speed PUBLIC
    Armourlib
    ${BOOST_LIBRARIES})

add_executable(BoxCollisionAvoidance_test 
    Constraints/TestBoxCollisionAvoidance.cpp)

target_link_libraries(BoxCollisionAvoidance_test PUBLIC
    trajlib
    IDlib
    Conslib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(BoxCollisionAvoidance_test PUBLIC
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

add_executable(MomentumRegressor_test 
    KinematicsDynamics/TestMomentumRegressor.cpp)

target_link_libraries(MomentumRegressor_test PUBLIC
    trajlib
    IDlib
    Conslib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(MomentumRegressor_test PUBLIC
    ${PINOCCHIO_FLAGS})

add_executable(InverseDynamics_test 
    KinematicsDynamics/TestInverseDynamics.cpp)

target_link_libraries(InverseDynamics_test PUBLIC
    trajlib
    IDlib
    Conslib
    ${BOOST_LIBRARIES} 
    pinocchio::pinocchio)

target_compile_options(InverseDynamics_test PUBLIC
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