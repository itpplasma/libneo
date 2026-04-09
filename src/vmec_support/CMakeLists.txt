add_library(vmec_support STATIC
    canonical_coordinates_mod.f90
    new_vmec_allocation_stuff.f90
    vmecinm_m.f90
    spline_vmec_data.f90
    vmec_field_tools.f90
)

target_link_libraries(vmec_support PUBLIC
    neo_util
    nctools
    interpolate
)

set_property(TARGET vmec_support PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::vmec_support ALIAS vmec_support)
