cmake_minimum_required(VERSION 3.10)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/field)

add_executable(test_biotsavart.x test_biotsavart.f90)
target_link_libraries(test_biotsavart.x PRIVATE neo
                                                util_for_test)
add_test(NAME test_biotsavart COMMAND test_biotsavart.x)

add_library(util_for_test_field OBJECT util_for_test_field.f90)
include_directories(${CMAKE_BINARY_DIR}/src/field)
target_link_libraries(util_for_test_field PRIVATE neo)

add_executable(test_field.x test_field.f90)
target_link_libraries(test_field.x PRIVATE neo
                                           util_for_test)
add_test(NAME test_field COMMAND test_field.x)

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

add_executable(test_polylag_field.x test_polylag_field.f90)
target_link_libraries(test_polylag_field.x PRIVATE neo
                                                  util_for_test_field
                                                  util_for_test)
add_test(NAME test_polylag_field COMMAND test_polylag_field.x)

set_tests_properties(
    test_biotsavart
    test_field
    test_example_field
    test_biotsavart_field
    test_polylag_field
    PROPERTIES LABELS "field")
