add_library(interpolate STATIC
    ../spl_three_to_five.f90
    acc_spline_build_order5.f90
    neo_bspline_base.f90
    neo_bspline_interp.f90
    neo_bspline_lsq.f90
    neo_bspline.f90
    batch_interpolate_types.f90
    batch_interpolate_1d.f90
    batch_interpolate_2d.f90
    batch_interpolate_3d.f90
    batch_interpolate.f90
    interpolate.f90)

if(LIBNEO_BUILD_TESTING)
    target_compile_definitions(interpolate PRIVATE LIBNEO_ENABLE_SPLINE_ORACLE=1)
endif()

# Apply OpenACC to interpolate library only (avoids nvfortran bugs in other modules)
# Link options must be INTERFACE so executables linking this library get -acc flag
if(ENABLE_OPENACC AND CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    # Memory model: "managed" keeps allocatable derived-type components such as
    # BatchSplineData%coeff coherent host<->device. Without it nvfortran does not
    # transfer the deep coeff array on "update self", so the host copy stays zero.
    # Must match the value used by consuming executables (e.g. SIMPLE).
    set(LIBNEO_OPENACC_MEM "managed" CACHE STRING
        "OpenACC memory model: managed|unified|separate")
    target_compile_options(interpolate PRIVATE -acc -gpu=ccnative,mem:${LIBNEO_OPENACC_MEM})
    target_link_options(interpolate INTERFACE -acc -gpu=ccnative,mem:${LIBNEO_OPENACC_MEM})
endif()

set_property(TARGET interpolate PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::interpolate ALIAS interpolate)
