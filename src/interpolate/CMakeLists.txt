add_library(interpolate OBJECT
    ../spl_three_to_five.f90
    batch_interpolate.f90
    interpolate.f90)

set_property(TARGET interpolate PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::interpolate ALIAS interpolate)
