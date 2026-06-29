set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/field)

add_executable(test_biotsavart.x test_biotsavart.f90)
target_link_libraries(test_biotsavart.x PRIVATE neo
                                                util_for_test)
add_test(NAME test_biotsavart COMMAND test_biotsavart.x)

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

add_executable(test_mgrid.x test_mgrid.f90)
target_link_libraries(test_mgrid.x PRIVATE neo
                                           util_for_test)
if(TARGET netcdf::netcdff)
    target_link_libraries(test_mgrid.x PRIVATE netcdf::netcdff)
elseif(TARGET netcdf::netcdf)
    target_link_libraries(test_mgrid.x PRIVATE netcdf::netcdf netcdff)
else()
    target_link_libraries(test_mgrid.x PRIVATE netcdf netcdff)
endif()
add_test(NAME test_mgrid COMMAND test_mgrid.x)

add_executable(test_stretch_coords.x test_stretch_coords.f90)
target_link_libraries(test_stretch_coords.x PRIVATE neo
                                                   magfie
                                                   util_for_test)
add_test(NAME test_stretch_coords COMMAND test_stretch_coords.x)


set_tests_properties(
    test_biotsavart
    test_example_field
    test_biotsavart_field
    test_mesh
    test_field_mesh
    test_polylag_field
    test_spline_field
    test_field
    test_mgrid
    test_stretch_coords
    PROPERTIES LABELS "field")
