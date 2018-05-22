file(GLOB SOURCE_FILES
        tree.c
        particle.c
        morton.c
        geometry.c)

add_executable(TreeNGB main.c ${SOURCE_FILES})

target_link_libraries(TreeNGB m)
