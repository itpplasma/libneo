set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/test/field)

add_executable(test_biotsavart.x test_biotsavart.f90)
target_link_libraries(test_biotsavart.x PRIVATE neo
                                                util_for_test)
add_test(NAME test_biotsavart COMMAND test_biotsavart.x)

add_executable(test_biotsavart_jacobian.x test_biotsavart_jacobian.f90)
target_link_libraries(test_biotsavart_jacobian.x PRIVATE neo)
add_test(NAME test_biotsavart_jacobian COMMAND test_biotsavart_jacobian.x)

add_executable(test_guiding_centre_field.x test_guiding_centre_field.f90)
target_link_libraries(test_guiding_centre_field.x PRIVATE neo)
add_test(NAME test_guiding_centre_field COMMAND test_guiding_centre_field.x)

add_executable(test_example_field.x test_example_field.f90)
target_link_libraries(test_example_field.x PRIVATE neo
                                                   util_for_test_field
                                                   util_for_test)
add_test(NAME test_example_field COMMAND test_example_field.x)

add_executable(test_biotsavart_field.x test_biotsavart_field.f90)
target_link_libraries(test_biotsavart_field.x PRIVATE neo
                                                      util_for_test_field
                                                      util_for_test)
add_test(NAME test_biotsavart_field COMMAND test_biotsavart_field.x)

add_executable(test_mesh.x test_mesh.f90)
target_link_libraries(test_mesh.x PRIVATE neo
                                          util_for_test)
add_test(NAME test_mesh COMMAND test_mesh.x)

add_executable(test_field_mesh.x test_field_mesh.f90)
target_link_libraries(test_field_mesh.x PRIVATE neo
                                               util_for_test)
add_test(NAME test_field_mesh COMMAND test_field_mesh.x)

add_executable(test_polylag_field.x test_polylag_field.f90)
target_link_libraries(test_polylag_field.x PRIVATE neo
                                                  util_for_test_field
                                                  util_for_test)
add_test(NAME test_polylag_field COMMAND test_polylag_field.x)

add_executable(test_spline_field.x test_spline_field.f90)
target_link_libraries(test_spline_field.x PRIVATE neo
                                              util_for_test_field
                                              util_for_test)
add_test(NAME test_spline_field COMMAND test_spline_field.x)

add_executable(test_field.x test_field.f90)
target_link_libraries(test_field.x PRIVATE neo
                                           util_for_test)
add_test(NAME test_field COMMAND test_field.x)

add_executable(test_fortfem_meridian.x test_fortfem_meridian.f90)
target_link_libraries(test_fortfem_meridian.x PRIVATE fortfem::api)
add_test(NAME test_fortfem_meridian COMMAND test_fortfem_meridian.x)

add_executable(test_free_boundary_coils.x test_free_boundary_coils.f90)
target_link_libraries(test_free_boundary_coils.x PRIVATE frozen_current fortio::fortio)
add_test(NAME test_free_boundary_coils COMMAND test_free_boundary_coils.x)

add_executable(test_frozen_current_field.x test_frozen_current_field.f90)
target_link_libraries(test_frozen_current_field.x PRIVATE frozen_current)
add_test(NAME test_frozen_current_field COMMAND test_frozen_current_field.x)
set_tests_properties(test_frozen_current_field PROPERTIES LABELS "field;slow")

add_executable(test_beams3d_field_grid.x test_beams3d_field_grid.f90)
target_link_libraries(test_beams3d_field_grid.x PRIVATE frozen_current neo)
add_test(NAME test_beams3d_field_grid COMMAND test_beams3d_field_grid.x)

add_executable(test_flux_label.x test_flux_label.f90)
target_link_libraries(test_flux_label.x PRIVATE frozen_current neo)
add_test(NAME test_flux_label COMMAND test_flux_label.x)

add_executable(test_frozen_current_limits.x test_frozen_current_limits.f90)
target_link_libraries(test_frozen_current_limits.x PRIVATE frozen_current)
add_test(NAME test_frozen_current_limits COMMAND test_frozen_current_limits.x)
set_tests_properties(test_frozen_current_limits PROPERTIES LABELS "field;slow")

add_executable(test_mgrid.x test_mgrid.f90)
target_link_libraries(test_mgrid.x PRIVATE neo
                                           util_for_test)
target_link_libraries(test_mgrid.x PRIVATE fortio::fortio)
add_test(NAME test_mgrid COMMAND test_mgrid.x)

add_executable(test_stretch_coords.x test_stretch_coords.f90)
target_link_libraries(test_stretch_coords.x PRIVATE neo
                                                   magfie
                                                   util_for_test)
add_test(NAME test_stretch_coords COMMAND test_stretch_coords.x)


set_tests_properties(
    test_biotsavart
    test_biotsavart_jacobian
    test_guiding_centre_field
    test_example_field
    test_biotsavart_field
    test_mesh
    test_field_mesh
    test_polylag_field
    test_spline_field
    test_field
    test_mgrid
    test_fortfem_meridian
    test_free_boundary_coils
    test_beams3d_field_grid
    test_flux_label
    test_stretch_coords
    PROPERTIES LABELS "field")

add_executable(test_cyl_grid_field.x test_cyl_grid_field.f90)
target_link_libraries(test_cyl_grid_field.x PRIVATE neo
                                                    util_for_test)
add_test(NAME test_cyl_grid_field COMMAND test_cyl_grid_field.x)

add_executable(test_makegrid_coils.x test_makegrid_coils.f90)
target_link_libraries(test_makegrid_coils.x PRIVATE neo
                                                    util_for_test)
add_test(NAME test_makegrid_coils COMMAND test_makegrid_coils.x)

add_executable(makegrid_coils_oracle.x makegrid_coils_oracle.f90)
target_link_libraries(makegrid_coils_oracle.x PRIVATE neo)
add_test(NAME test_makegrid_coils_python_oracle
         COMMAND ${Python_EXECUTABLE}
                 ${PROJECT_SOURCE_DIR}/test/scripts/compare_makegrid_coils.py
                 --exe $<TARGET_FILE:makegrid_coils_oracle.x>
                 --work-dir ${PROJECT_BINARY_DIR}/test/field/makegrid_oracle)
set_tests_properties(test_makegrid_coils test_makegrid_coils_python_oracle
                     PROPERTIES LABELS "field")
