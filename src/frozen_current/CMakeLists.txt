# Frozen-current field: equilibrium sources, free-boundary coils and the
# per-mode FortFEM solve. Separate from neo so that only its consumers build
# FortFEM.
add_library(frozen_current STATIC
    frozen_current_leaves.f90
    equilibrium_source.f90
    vmec_equilibrium_source.f90
    free_boundary_coils.f90
    frozen_current_field.f90
)
target_link_libraries(frozen_current PUBLIC vmec_support neo_field fortnum fortfem::api)
target_link_libraries(frozen_current PRIVATE fortio::fortio)
if(ENABLE_OPENMP)
    target_link_libraries(frozen_current PUBLIC OpenMP::OpenMP_Fortran)
endif()
set_property(TARGET frozen_current PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::frozen_current ALIAS frozen_current)
