cmake_minimum_required(VERSION 3.10)

add_library(neo_field OBJECT
    field_base.f90
    example_field.f90
    biotsavart.f90
    biotsavart_field.f90
    field_mesh.f90
    ../polylag/polylag_5.f90
    polylag_field.f90
    field.f90)

# target_link_libraries(neo_field PRIVATE $<TARGET_OBJECTS:neo_polylag>)
# target_include_directories(neo_field PRIVATE ${CMAKE_BINARY_DIR}/src/polylag)
