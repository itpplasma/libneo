set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/field)

add_library(util_for_test_field OBJECT util_for_test_field.f90)
target_link_libraries(util_for_test_field PRIVATE neo)

add_library(util_for_test_jorek_field OBJECT util_for_test_jorek_field.f90)


add_executable(setup_test_jorek_field.x setup_test_jorek_field.f90)
target_link_libraries(setup_test_jorek_field.x PRIVATE neo
                                                       hdf5_tools
                                                       util_for_test
                                                       util_for_test_jorek_field)
add_test(NAME setup_test_jorek_field COMMAND setup_test_jorek_field.x)
set_tests_properties(setup_test_jorek_field PROPERTIES 
                     FIXTURES_SETUP "SetupTestFieldJorek")

add_executable(cleanup_test_jorek_field.x cleanup_test_jorek_field.f90)
target_link_libraries(cleanup_test_jorek_field.x PRIVATE neo
                                                         hdf5_tools
                                                         util_for_test
                                                         util_for_test_jorek_field)
add_test(NAME cleanup_test_jorek_field COMMAND cleanup_test_jorek_field.x)
set_tests_properties(cleanup_test_jorek_field PROPERTIES 
                     FIXTURES_CLEANUP "SetupTestFieldJorek")


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

add_executable(test_jorek_field.x test_jorek_field.f90)
target_link_libraries(test_jorek_field.x PRIVATE neo
                                                 hdf5_tools
                                                 util_for_test
                                                 util_for_test_jorek_field)
add_test(NAME test_jorek_field COMMAND test_jorek_field.x)
set_tests_properties(test_jorek_field PROPERTIES 
                     FIXTURES_REQUIRED "SetupTestFieldJorek")

add_executable(test_field.x test_field.f90)
target_link_libraries(test_field.x PRIVATE neo
                                           util_for_test
                                           util_for_test_jorek_field)
add_test(NAME test_field COMMAND test_field.x)
set_tests_properties(test_field PROPERTIES 
                     FIXTURES_REQUIRED "SetupTestFieldJorek")


set_tests_properties(
    test_biotsavart
    test_example_field
    test_biotsavart_field
    test_mesh
    test_field_mesh
    test_polylag_field
    test_spline_field
    test_jorek_field
    test_field
    PROPERTIES LABELS "field")
