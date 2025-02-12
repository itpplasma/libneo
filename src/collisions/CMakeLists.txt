include(FetchContent)

# Try to find fortran_stdlib first
find_package(fortran_stdlib REQUIRED)

if(NOT fortran_stdlib_FOUND)
    message(STATUS "fortran-stdlib not found, fetching from GitHub...")

    FetchContent_Declare(
        fortran_stdlib
        GIT_REPOSITORY https://github.com/fortran-lang/stdlib.git
        GIT_TAG main  # You may want to specify a stable version
    )

    FetchContent_MakeAvailable(fortran_stdlib)
endif()

## Add the library to your project
#add_subdirectory(${fortran_stdlib_SOURCE_DIR})
#add_subdirectory(${fortran_stdlib_BINARY_DIR})

add_library(collision_freqs STATIC
    collision_freqs.f90)

set_property(TARGET collision_freqs PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::collision_freqs ALIAS collision_freqs)
target_link_libraries(collision_freqs PUBLIC fortran_stdlib::fortran_stdlib)

