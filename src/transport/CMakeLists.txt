add_library(transport STATIC
    transport.f90
    gen_laguerre_rule.f90)

set_property(TARGET transport PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})

add_library(LIBNEO::transport ALIAS transport)

add_dependencies(transport collision_freqs)

target_link_libraries(transport PUBLIC collision_freqs)

