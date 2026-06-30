add_library(neo_util STATIC
    libneo_kinds.f90
    math_constants.f90
    util.f90
    binsrc.f90
    plag_coeff.f90
    arnoldi.f90
    rusage_type.f90
    system_utility.f90
    simpson_integration.f90
    local_rusage.c
)

set_property(TARGET neo_util PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::util ALIAS neo_util)
