file(GLOB SOURCE_FILES
        tree.c
        main.c
        particle.c
        morton.c)

add_executable(TreeNGB ${SOURCE_FILES})

target_link_libraries(TreeNGB m)
