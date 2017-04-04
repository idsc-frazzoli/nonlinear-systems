cmake_minimum_required(VERSION 2.6)
project(dynamic_system_simulator)

add_executable(dynamic_system_simulator main.cpp)

install(TARGETS dynamic_system_simulator RUNTIME DESTINATION bin)
