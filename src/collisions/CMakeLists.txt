find_package(fortran_stdlib QUIET)

if (NOT fortran_stdlib_FOUND)
    message(STATUS "fortran_stdlib not found, building it from source")
    include(ExternalProject)

    set(prefix "${CMAKE_BINARY_DIR}/fortran_stdlib")
    set(FORTRAN_STDLIB_LIBRARY ${prefix}/src/fortran_stdlib-build/src/libfortran_stdlib.a)
    set(FORTRAN_STDLIB_INCLUDE_DIR ${prefix}/src/fortran_stdlib-build/src/mod_files)

    execute_process(
    COMMAND ${CMAKE_COMMAND} -E make_directory ${FORTRAN_STDLIB_INCLUDE_DIR}
    )

    ExternalProject_Add(
        fortran_stdlib
        GIT_REPOSITORY https://github.com/fortran-lang/stdlib.git
        PREFIX ${prefix}
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/fortran_stdlib
        CMAKE_ARGS -DBUILD_TESTING="off"
        BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> #--target fortran_stdlib/install
        INSTALL_COMMAND ""
        BUILD_BYPRODUCTS ${FORTRAN_STDLIB_LIBRARY}
            ${FORTRAN_STDLIB_INCLUDE_DIR}
    )

    add_library(fortran_stdlib::fortran_stdlib STATIC IMPORTED)
    set_target_properties(fortran_stdlib::fortran_stdlib PROPERTIES
        IMPORTED_LOCATION ${FORTRAN_STDLIB_LIBRARY}
        INTERFACE_INCLUDE_DIRECTORIES ${FORTRAN_STDLIB_INCLUDE_DIR}
    )
endif()

add_library(collision_freqs STATIC
    collision_freqs.f90)

set_property(TARGET collision_freqs PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::collision_freqs ALIAS collision_freqs)
target_link_libraries(collision_freqs PUBLIC fortran_stdlib::fortran_stdlib)
target_include_directories(collision_freqs PUBLIC ${FORTRAN_STDLIB_INCLUDE_DIR})