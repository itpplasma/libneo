add_library(interpolate OBJECT
    interpolate.f90)

set_property(TARGET interpolate PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::interpolate ALIAS interpolate)
