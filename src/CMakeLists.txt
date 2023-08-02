# This file provides the list of source files for libneo.
# The file CMakeLists.txt contains the main configuration code.
# Main goal (of the created build file): creating the library file.
add_library(magfie SHARED
    boozer.f90
    io.f90
    libneo_kinds.f90
    math_constants.f90

    # VMEC
    canonical_coordinates_mod.f90
    nctools_module.f90
    new_vmec_allocation_stuff.f90
    spline_vmec_data.f90
    vmecinm_m.f90
    spl_three_to_five.f90
    plag_coeff.f90
    odeint_rkf45.f90
    binsrc.f90
    get_canonical_coordinates.f90
    boozer_converter.f90
    magfie/magfie_vmec.f90

    # BDIVFREE
    magfie/amn_mod.f90
    magfie/bdivfree.f90
    magfie/bdivfree_mod.f90
    magfie/inthecore_mod.f90
    magfie/theta_rz_mod.f90
    magfie/input_files.f90
    magfie/field_eq_mod.f90
    magfie/extract_fluxcoord_mod.f90
    magfie/field_divB0.f90
    magfie/period_mod.f90
    magfie/field_c_mod.f90
    magfie/field_mod.f90
    magfie/spline5_RZ.f90
    magfie/field_mesh3d.f90
    magfie/magfield_mod.f90
    magfie/polylag_3.f
    magfie/magfie_cyl.f90
)
set_source_files_properties(boozer.f90 PROPERTIES COMPILE_OPTIONS "-Wno-error=unused-variable;-Wno-error=maybe-uninitialized")
set_source_files_properties(field_divB0.f90 PROPERTIES COMPILE_OPTIONS "-Wno-error=compare-reals;-Wno-error=conversion")
set_source_files_properties(field_mesh3d.f90 PROPERTIES COMPILE_OPTIONS "-Wno-error=compare-reals")

add_executable(standardise_geqdsk.x
    magfie/standardise_geqdsk.f90
    magfie/geqdsk_tools.f90
    plag_coeff.f90
    libneo_kinds.f90
    math_constants.f90
)

target_link_libraries(magfie
  rkf45 netcdf netcdff ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
