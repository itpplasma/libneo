set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/spectre)

add_executable(test_spectre_reader.x test_spectre_reader.f90)
target_link_libraries(test_spectre_reader.x PRIVATE spectre_support
                                                    hdf5_tools
                                                    util_for_test)
add_test(NAME test_spectre_reader
         COMMAND test_spectre_reader.x
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/spectre_test.h5)

add_executable(test_spectre_basis.x test_spectre_basis.f90)
target_link_libraries(test_spectre_basis.x PRIVATE spectre_support
                                                   hdf5_tools
                                                   util_for_test)
add_test(NAME test_spectre_basis
         COMMAND test_spectre_basis.x
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/spectre_test.h5
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/reference_eval.dat)
