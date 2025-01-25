add_library(neo_field STATIC
    field_base.f90
    example_field.f90
    circular_tokamak_field.f90
    biotsavart.f90
    biotsavart_field.f90
    field_mesh.f90
    polylag_field.f90
    spline_field.f90
    jorek_field.f90
    mesh.f90
    mesh_2d.f90
    field.f90)

set_property(TARGET neo_field PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::field ALIAS neo_field)
