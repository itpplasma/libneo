cmake_minimum_required(VERSION 3.10)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/polylag)

add_executable(test_polylag_5.x test_polylag_5.f90)
target_link_libraries(test_polylag_5.x PRIVATE neo
                                               # neo_polylag
                                               util_for_test)
target_include_directories(test_polylag_5.x PRIVATE ${CMAKE_BINARY_DIR}/src/field)
add_test(NAME test_polylag_5 COMMAND test_polylag_5.x)

set_tests_properties(
    test_polylag_5
    PROPERTIES LABELS "polylag")
