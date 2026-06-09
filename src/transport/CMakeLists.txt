add_library(transport STATIC
    transport.f90)

set_property(TARGET transport PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})

add_library(LIBNEO::transport ALIAS transport)
target_link_libraries(transport PUBLIC collision_freqs species neo_math)

