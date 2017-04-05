cmake_minimum_required(VERSION 2.8)
project(dynamic_system_simulator)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -std=c++0x -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -std=c++0x")#for debug

include_directories(${PROJECT_SOURCE_DIR}/include/spline)
include_directories(${PROJECT_SOURCE_DIR}/include/numerical_integration)

add_executable(dynamic_system_simulator main.cpp)
install(TARGETS dynamic_system_simulator RUNTIME DESTINATION bin)

