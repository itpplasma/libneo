add_library(neo_polylag OBJECT
    polylag_5.f90)

set_property(TARGET neo_polylag PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::polylag ALIAS neo_polylag)
