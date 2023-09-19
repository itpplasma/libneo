# This file provides the list of source files for libneo.
# The file CMakeLists.txt contains the main configuration code.
# Main goal (of the created build file): creating the library file.

add_library(mc_efit SHARED
    libneo_kinds.f90
    math_constants.f90
    spl_three_to_five.f90
    magfie/polylag_3.f
    magfie/polylag_5.f
    magfie/spline5_RZ.f90
    magfie/amn_mod.f90
    magfie/bdivfree.f90
    magfie/bdivfree_mod.f90
    magfie/extract_fluxcoord_mod.f90
    magfie/field_c_mod.f90
    magfie/field_divB0.f90
    magfie/field_eq_mod.f90
    magfie/field_mesh3d.f90
    magfie/field_mod.f90
    magfie/input_files.f90
    magfie/inthecore_mod.f90
    magfie/magfie_cyl.f90
    magfie/magfield_mod.f90
    magfie/neo_magfie_mod.f90
    magfie/period_mod.f90
    magfie/theta_rz_mod.f90
    contrib/rkf45.f90
    odeint_rkf45.f90
    MC/alpha_lifetime_mod.f90
    MC/chamb_divB0.f90
    MC/rng.f90
    MC/collis_nbi.f90
    MC/sub_alpha_lifetime.f90
)


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
    magfie/geqdsk_tools.f90
    magfie/coil_tools.f90
    fftw3.f90
)
set_source_files_properties(boozer.f90 PROPERTIES COMPILE_OPTIONS "-Wno-error=unused-variable;-Wno-error=maybe-uninitialized")
set_source_files_properties(field_divB0.f90 PROPERTIES COMPILE_OPTIONS "-Wno-error=compare-reals;-Wno-error=conversion")
set_source_files_properties(field_mesh3d.f90 PROPERTIES COMPILE_OPTIONS "-Wno-error=compare-reals")

add_executable(standardise_geqdsk.x
    magfie/standardise_geqdsk.f90
    magfie/geqdsk_tools.f90
    plag_coeff.f90
    math_constants.f90
    libneo_kinds.f90
)
add_executable(vacfield.x
    magfie/vacfield.f90
    magfie/coil_tools.f90
    fftw3.f90
    math_constants.f90
    libneo_kinds.f90
)

find_library(fftw_lib NAMES "fftw3")
find_path(fftw_include_dir NAMES "fftw3.h")

find_package(HDF5 COMPONENTS C Fortran Fortran_HL REQUIRED)
if (NOT ${HDF5_FOUND})
	message(FATAL_ERROR "HDF5 not found, please set environment variable HDF5_ROOT to the correct path")
endif ()
include_directories(${HDF5_INCLUDE_DIRS})

set(hdf5_tools_include_dir ${CMAKE_BINARY_DIR})

target_link_libraries(magfie
  rkf45
  netcdf
  netcdff
  hdf5_tools
  ${fftw_lib}
  ${HDF5_HL_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
)

target_link_libraries(vacfield.x
  hdf5_tools
  ${HDF5_HL_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${fftw_lib}
)

target_include_directories(magfie PUBLIC ${hdf5_tools_include_dir} ${fftw_include_dir})
target_include_directories(vacfield.x PUBLIC ${hdf5_tools_include_dir} ${fftw_include_dir})
