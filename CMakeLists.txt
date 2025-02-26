# Minimum required CMake version
cmake_minimum_required(VERSION 3.18)

# Project name
project(RAPTOR
    DESCRIPTION "Rapid and Robust Trajectory Optimization for Robots")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include(CTest)
set(CTEST_OUTPUT_ON_FAILURE TRUE)
set(CTEST_CUSTOM_TESTS_IGNORE "")

set(CTEST_JUNIT_OUTPUT_PATH "${CMAKE_BINARY_DIR}/TestResults")
set(CTEST_CUSTOM_TESTS_IGNORE "")

enable_testing()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp")

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Number of threads that openmp will use
if(DEFINED NUM_THREADS)
  add_definitions(-DNUM_THREADS=${NUM_THREADS})
  message(STATUS "Number of threads used in openmp: ${NUM_THREADS}")
else()
    # Get the maximum number of threads using nproc
    execute_process(
        COMMAND nproc
        OUTPUT_VARIABLE MAX_THREADS
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    add_definitions(-DNUM_THREADS=${MAX_THREADS})
    message(STATUS "Maximize number of threads used in openmp: ${MAX_THREADS}")
endif()

# Find Python binding library
find_package(Python 3.10
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# Find Boost library
find_package(Boost REQUIRED COMPONENTS 
    system
    filesystem
    serialization
)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

# Find GSL library
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})

# Find Eigen library
find_package(Eigen3 3.3.7 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

# Find Pinocchio library
find_package(pinocchio REQUIRED)
include_directories(${PINOCCHIO_INCLUDE_DIRS})
link_directories(${PINOCCHIO_LIBRARY_DIRS})

# Find ompl library
find_package(ompl REQUIRED)
set(OMPL_INCLUDE_DIRS /usr/local/include/ompl-1.6)
set(OMPL_LIBRARY_DIRS /usr/local/share)
include_directories(${OMPL_INCLUDE_DIRS})
link_directories(${OMPL_LIBRARY_DIRS})

set(PINOCCHIO_FLAGS
    ${pinocchio_CFLAGS_OTHER}
    -Wno-ignored-attributes
    -Wno-invalid-partial-specialization    # to silence warning with unsupported Eigen Tensor
    -DPINOCCHIO_URDFDOM_TYPEDEF_SHARED_PTR
    -DPINOCCHIO_URDFDOM_USE_STD_SHARED_PTR
)
# Find YAML library
find_package(yaml-cpp REQUIRED)
include_directories(${YAML_CPP_INCLUDE_DIR})
link_directories(${YAML_CPP_LIBRARY_DIR})

# Find Ipopt library
set(IPOPT_INCLUDE_DIR /usr/local/include/coin-or)
set(IPOPT_LIBRARY_DIR /usr/local/lib)
include_directories(${IPOPT_INCLUDE_DIR})
link_directories(${IPOPT_LIBRARY_DIR})

# Find HSL library
set(HSL_INCLUDE_DIR /usr/local/include/coin/ThirdParty)
set(HSL_LIBRARY_DIR /usr/local/lib)
include_directories(${HSL_INCLUDE_DIR})
link_directories(${HSL_LIBRARY_DIR})

# # Find pardiso library
# set(PARDISO_INCLUDE_DIR /usr/local/include)
# set(PARDISO_LIBRARY_DIR /usr/local/lib)
# include_directories(${PARDISO_INCLUDE_DIR})
# link_directories(${PARDISO_LIBRARY_DIR})
# include_directories(Utils/include)

# Find MOSEK library
# set(MOSEK_INCLUDE_DIR /opt/mosek/10.2/tools/platform/linux64x86/h)
# set(MOSEK_LIBRARY_DIR /opt/mosek/10.2/tools/platform/linux64x86/bin)
# include_directories(${MOSEK_INCLUDE_DIR})
# link_directories(${MOSEK_LIBRARY_DIR})

add_subdirectory(Trajectories)
add_subdirectory(KinematicsDynamics)
add_subdirectory(Optimization)
add_subdirectory(Costs)
add_subdirectory(Constraints)
add_subdirectory(Examples)
add_subdirectory(Tests)

# include
include_directories(
    ${PROJECT_SOURCE_DIR}/Costs/include
    ${PROJECT_SOURCE_DIR}/Constraints/include
    ${PROJECT_SOURCE_DIR}/KinematicsDynamics/include
    ${PROJECT_SOURCE_DIR}/Optimization/include
    ${PROJECT_SOURCE_DIR}/Trajectories/include
    ${PROJECT_SOURCE_DIR}/Utils/include
)

## Install
install(TARGETS 
    trajlib
    IDlib
    Costslib
    Conslib
    Optlib
    Kinovalib
    Digitlib
    Taloslib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)