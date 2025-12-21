add_library(interpolate STATIC
    ../spl_three_to_five.f90
    acc_tridiag_pcr.f90
    acc_spline_build_order5.f90
    batch_interpolate_types.f90
    batch_interpolate_1d.f90
    batch_interpolate_2d.f90
    batch_interpolate_3d.f90
    batch_interpolate.f90
    interpolate.f90)

set_property(TARGET interpolate PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::interpolate ALIAS interpolate)
