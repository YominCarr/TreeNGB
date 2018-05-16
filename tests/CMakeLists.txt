include_directories (
        ${CMAKE_INCLUDE_PATH}
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/external/gmock-1.7.0/
        ${CMAKE_SOURCE_DIR}/external/gtest-1.7.0/
)

file(GLOB SOURCE_FILES
        ../src/tree.c
        ../src/particle.c
        ../src/morton.c)

file(GLOB TEST_SOURCES
        testMorton.cpp
        testTree.cpp
        testParticle.cpp
)
list(APPEND TEST_SOURCES "allTests.cpp")
add_executable (TreeNGB-tests ${SOURCE_FILES} ${TEST_SOURCES})
target_link_libraries (TreeNGB-tests gmock gtest)

message(STATUS ${SOURCE_FILES})
