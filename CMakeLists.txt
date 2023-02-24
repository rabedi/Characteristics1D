cmake_minimum_required(VERSION 3.10)
project(MOC1D VERSION 0.1.0
        DESCRIPTION "One-dimensional Method of Characteristics Solver"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

aux_source_directory(. DIR_UTIL_SRCS)

add_executable(solver main.cpp ${DIR_UTIL_SRCS})


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

