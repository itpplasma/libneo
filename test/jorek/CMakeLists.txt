set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/jorek)

add_executable(test_jorek_restart.x test_jorek_restart.f90)
target_link_libraries(test_jorek_restart.x PRIVATE jorek_support
                                                   hdf5_tools
                                                   util_for_test)
add_test(NAME test_jorek_restart COMMAND test_jorek_restart.x)

add_executable(test_jorek_restart_fixture.x test_jorek_restart_fixture.f90)
target_link_libraries(test_jorek_restart_fixture.x PRIVATE jorek_support
                                                           hdf5_tools
                                                           util_for_test)
add_test(NAME test_jorek_restart_fixture COMMAND test_jorek_restart_fixture.x)
