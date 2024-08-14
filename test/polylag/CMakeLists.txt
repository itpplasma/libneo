set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/polylag)

add_executable(test_polylag_5.x test_polylag_5.f90)
target_link_libraries(test_polylag_5.x PRIVATE neo
                                               util_for_test)
add_test(NAME test_polylag_5 COMMAND test_polylag_5.x)

set_tests_properties(
    test_polylag_5
    PROPERTIES LABELS "polylag")
