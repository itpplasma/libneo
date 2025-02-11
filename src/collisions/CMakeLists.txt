add_library(collisions STATIC
    perp_coll_frequency.f90)

set_property(TARGET collisions PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::collisions ALIAS collisions)

