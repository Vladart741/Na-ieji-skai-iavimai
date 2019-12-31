cmake_minimum_required(VERSION 3.14)
project(NS)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
set(SOURCE_FILES main.cpp)
add_executable(NS main.cpp)