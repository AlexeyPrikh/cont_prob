cmake_minimum_required(VERSION 3.17)
project(heat_3D)

set(CMAKE_CXX_STANDARD 17)

# added -fopenmp
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")

#add_executable(heat_3D main.cpp Solution.cpp Solution.h Special_functions.h Special_functions.cpp solver_3D.cpp solver_3D.h)

add_executable(heat_3D main.cpp Solution.cpp Special_functions.cpp solver_3D.cpp)