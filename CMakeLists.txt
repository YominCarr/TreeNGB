cmake_minimum_required(VERSION 2.8)
project(TreeNGB)

option (BUILD_TESTS "Build tests" ON)

if(BUILD_TESTS)
    add_definitions(-DBUILD_TESTS)
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O0 -Wall -fopenmp")

message(STATUS "Compiling with ${CMAKE_C_FLAGS}")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

add_subdirectory (src)
add_subdirectory (external)
if(BUILD_TESTS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -Wall -fopenmp")
    add_subdirectory (tests)
endif()
