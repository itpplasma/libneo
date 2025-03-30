add_library(collision_freqs STATIC
    collision_freqs.f90)

find_package(GSL REQUIRED)

set_property(TARGET collision_freqs PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})

add_library(LIBNEO::collision_freqs ALIAS collision_freqs)
add_dependencies(collision_freqs species)

target_link_libraries(collision_freqs PUBLIC GSL::gsl)
target_link_libraries(collision_freqs PUBLIC neo
                                            species)
target_link_directories(collision_freqs PUBLIC species)
