cmake_minimum_required(VERSION 2.8)
project(TreeNGB)

set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -std=c99 -O2 -Wall -fopenmp")


message(STATUS "Compiling with ${CMAKE_C_FLAGS}")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

set(SOURCE_FILES
        src/tree.c
        src/main.c)

add_executable(TreeNGB ${SOURCE_FILES})
