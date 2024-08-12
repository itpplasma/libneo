cmake_minimum_required(VERSION 3.10)

add_library(util_for_test OBJECT util_for_test.f90)

add_executable(test_biotsavart.x test_biotsavart.f90)
target_link_libraries(test_biotsavart.x PRIVATE libneo_field
                                                util_for_test)

add_library(util_for_test_field OBJECT util_for_test_field.f90)
include_directories(${CMAKE_BINARY_DIR}/src/field)
target_link_libraries(util_for_test_field PRIVATE libneo_field)

add_executable(test_field.x test_field.f90)
target_link_libraries(test_field.x PRIVATE libneo_field
                                           util_for_test)

add_executable(test_example_field.x test_example_field.f90)
target_link_libraries(test_example_field.x PRIVATE libneo_field
                                                   util_for_test_field
                                                   util_for_test)

add_executable(test_biotsavart_field.x test_biotsavart_field.f90)
target_link_libraries(test_biotsavart_field.x PRIVATE libneo_field
                                                      util_for_test_field
                                                      util_for_test)
