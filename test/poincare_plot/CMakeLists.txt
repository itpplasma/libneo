set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/poincare_plot)

add_executable(test_poincare_plot.x test_poincare_plot.f90)
target_link_libraries(test_poincare_plot.x PRIVATE neo
                                               util_for_test)
add_test(NAME test_poincare_plot COMMAND test_poincare_plot.x)

set_tests_properties(
    test_poincare_plot
    PROPERTIES LABELS "poincare_plot")
