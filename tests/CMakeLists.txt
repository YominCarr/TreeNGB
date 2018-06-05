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
        ../src/morton.c
        ../src/geometry.c)

file(GLOB TEST_SOURCES
        testMorton.cpp
        testTree.cpp
        testParticle.cpp
        testGeometry.cpp
)
list(APPEND TEST_SOURCES "allTests.cpp")
add_executable (TreeNGB-tests ${SOURCE_FILES} ${TEST_SOURCES})
target_link_libraries (TreeNGB-tests gmock gtest)

add_custom_target(TreeNGB-profile COMMAND pwd COMMAND valgrind --tool=callgrind ./TreeNGB-tests --gtest_filter=*complexNeighbourFindingAgainstBruteForce* WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}" DEPENDS TreeNGB-tests)

message(STATUS ${SOURCE_FILES})
