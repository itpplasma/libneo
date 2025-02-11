add_library(collision_freqs STATIC
    collision_freqs.f90)

set_property(TARGET collision_freqs PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::collision_freqs ALIAS collision_freqs)

