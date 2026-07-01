# field_base: abstract field interface (neo_field_base). No internal libneo
# deps, so it sits at Tier 0 with the other leaves.
add_library(field_base STATIC
    field_base.f90)
set_property(TARGET field_base PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::field_base ALIAS field_base)

# neo_field: concrete field implementations built on field_base, plus the mesh
# containers they consume.
add_library(neo_field STATIC
    mesh.f90
    mesh_2d.f90
    example_field.f90
    circular_tokamak_field.f90
    biotsavart.f90
    biotsavart_field.f90
    field_mesh.f90
    polylag_field.f90
    spline_field.f90
    mgrid.f90
    field.f90)

target_link_libraries(neo_field PUBLIC field_base)
target_link_libraries(neo_field PRIVATE
    interpolate
    neo_polylag
)
if(TARGET netcdf::netcdff)
    target_link_libraries(neo_field PRIVATE netcdf::netcdff)
elseif(TARGET netcdf::netcdf)
    target_link_libraries(neo_field PRIVATE netcdf::netcdf netcdff)
else()
    target_link_libraries(neo_field PRIVATE netcdf netcdff)
endif()
set_property(TARGET neo_field PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::field ALIAS neo_field)
