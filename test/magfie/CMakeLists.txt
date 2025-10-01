set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/magfie)

# Test executable for Biot-Savart Fourier routines
add_executable(test_coil_tools_biot_savart.x test_coil_tools_biot_savart.f90)
target_link_libraries(test_coil_tools_biot_savart.x PRIVATE 
    neo
    magfie
    util_for_test
    ${FFTW3_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
)

# Register the test
add_test(NAME test_coil_tools_biot_savart COMMAND test_coil_tools_biot_savart.x)

# Set test properties
set_tests_properties(test_coil_tools_biot_savart PROPERTIES
    LABELS "magfie;biot_savart"
    TIMEOUT 300
)

# Integration test: Compare GPEC Fourier vs coil_tools vector_potential on AUG data
add_test(NAME test_aug_biot_savart_comparison
    COMMAND ${CMAKE_COMMAND}
        -DVACFIELD_EXE=$<TARGET_FILE:vacfield.x>
        -DTEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/test_data
        -DPROJECT_SOURCE_DIR=${CMAKE_SOURCE_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/run_aug_comparison.cmake
)

set_tests_properties(test_aug_biot_savart_comparison PROPERTIES
    LABELS "magfie;biot_savart;integration"
    TIMEOUT 600
    DEPENDS vacfield.x
)

add_executable(benchmark_biot_savart.x benchmark_biot_savart.f90)
target_link_libraries(benchmark_biot_savart.x PRIVATE
    neo
    magfie
    ${FFTW3_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
)

add_test(NAME benchmark_biot_savart
    COMMAND benchmark_biot_savart.x ${CMAKE_SOURCE_DIR}/test/magfie/test_data
)

set_tests_properties(benchmark_biot_savart PROPERTIES
    LABELS "magfie;biot_savart;benchmark"
)

# Integration test: single-coil scenario with analytical verification
add_test(NAME test_single_coil_biot_savart
    COMMAND ${CMAKE_COMMAND}
        -DVACFIELD_EXE=$<TARGET_FILE:vacfield.x>
        -DTEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/test_data
        -DPROJECT_SOURCE_DIR=${CMAKE_SOURCE_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/run_single_coil_case.cmake
)

set_tests_properties(test_single_coil_biot_savart PROPERTIES
    LABELS "magfie;biot_savart;integration;single_coil"
    TIMEOUT 600
    DEPENDS vacfield.x
)
