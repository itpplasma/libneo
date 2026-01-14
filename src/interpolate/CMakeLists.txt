add_library(interpolate STATIC
    ../spl_three_to_five.f90
    acc_spline_build_order5.f90
    batch_interpolate_types.f90
    batch_interpolate_1d.f90
    batch_interpolate_2d.f90
    batch_interpolate_3d.f90
    batch_interpolate.f90
    cgls_small.f90
    cgls_dense.f90
    interpolate.f90)

if(LIBNEO_BUILD_TESTING)
    target_compile_definitions(interpolate PRIVATE LIBNEO_ENABLE_SPLINE_ORACLE=1)
endif()

# Apply OpenACC to interpolate library only (avoids nvfortran bugs in other modules)
# Link options must be INTERFACE so executables linking this library get -acc flag
if(ENABLE_OPENACC AND CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    target_compile_options(interpolate PRIVATE -acc -gpu=ccnative)
    target_link_options(interpolate INTERFACE -acc -gpu=ccnative)
endif()

set_property(TARGET interpolate PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::interpolate ALIAS interpolate)
