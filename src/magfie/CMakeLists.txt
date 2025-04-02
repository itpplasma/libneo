add_library(magfie STATIC
    # BDIVFREE
    amn_mod.f90
    bdivfree.f90
    bdivfree_mod.f90
    inthecore_mod.f90
    theta_rz_mod.f90
    input_files.f90
    field_eq_mod.f90
    extract_fluxcoord_mod.f90
    field_divB0.f90
    period_mod.f90
    field_c_mod.f90
    field_mod.f90
    spline5_RZ.f90
    field_mesh3d.f90
    magfield_mod.f90
    polylag_3.f
    magfie_cyl.f90
    geqdsk_tools.f90
    coil_tools.f90
    ../fftw3.F90
)
target_link_libraries(magfie PRIVATE
    neo
    rkf45
    netcdf
    netcdff
    hdf5_tools
    FFTW::DoubleThreads
    FFTW::Double
)
set_property(TARGET neo_field PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::magfie ALIAS magfie)

