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