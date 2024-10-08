set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/poincare)

add_executable(test_poincare.x test_poincare.f90)
target_link_libraries(test_poincare.x PRIVATE neo
                                               util_for_test
                                               util_for_test_jorek_field)
add_test(NAME test_poincare COMMAND test_poincare.x)
set_tests_properties(test_poincare PROPERTIES 
                     FIXTURES_REQUIRED "SetupTestFieldJorek")

set_tests_properties(
    test_poincare
    PROPERTIES LABELS "poincare")
