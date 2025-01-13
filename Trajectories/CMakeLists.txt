add_library(trajlib SHARED
    src/Trajectories.cpp
    src/FourierCurves.cpp 
    src/FixedFrequencyFourierCurves.cpp
    src/BezierCurves.cpp 
    src/ArmourBezierCurves.cpp
    src/PiecewiseBezierCurves.cpp
    src/Polynomials.cpp 
    src/Plain.cpp
    src/TrajectoryData.cpp
    src/TrajectoryGroup.cpp
)

target_include_directories(trajlib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/Utils/include
)