# Minimum required CMake version
cmake_minimum_required(VERSION 3.12)

# Project name
project(IDTO
    DESCRIPTION "Inverse Dynamics Trajectory Optimization")

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)

# Find Eigen library
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

# Find GSL library
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

## Trajectories library
add_library(trajlib SHARED
    Trajectories/src/FourierCurves.cpp 
    Trajectories/src/Trajectories.cpp)

target_include_directories(trajlib PUBLIC
    Trajectories/include)

## Example scripts
add_executable(IDTO_example 
    Examples/Digit/DigitSingleStep.cpp)

target_include_directories(IDTO_example PRIVATE
    Trajectories/include)

target_link_libraries(IDTO_example PRIVATE
    trajlib 
    ${EIGEN3_LIBRARIES} 
    ${GSL_LIBRARIES})

## Install
install(TARGETS trajlib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(TARGETS IDTO_example
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
