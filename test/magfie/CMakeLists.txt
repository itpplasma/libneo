set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/magfie)

find_package(Python REQUIRED COMPONENTS Interpreter)

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

add_executable(test_coil_tools_vector_potential_derivs.x
    test_coil_tools_vector_potential_derivs.f90)
target_link_libraries(test_coil_tools_vector_potential_derivs.x PRIVATE
    neo
    magfie
    util_for_test
    ${FFTW3_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
)
target_compile_definitions(test_coil_tools_vector_potential_derivs.x PRIVATE
    COIL_TOOLS_PLOT_DIR="${CMAKE_BINARY_DIR}/test/magfie")

add_test(NAME test_coil_tools_vector_potential_derivs
    COMMAND test_coil_tools_vector_potential_derivs.x)
set_tests_properties(test_coil_tools_vector_potential_derivs PROPERTIES
    LABELS "magfie;coil_tools;derivatives"
    TIMEOUT 300
)

add_test(NAME plot_vector_potential_deriv_errors
    COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/python/scripts/plot_derivative_error.py
        --r-input ${CMAKE_BINARY_DIR}/test/magfie/coil_tools_dAphi_dR_error.dat
        --r-output ${CMAKE_BINARY_DIR}/test/magfie/coil_tools_dAphi_dR_error.png
        --z-input ${CMAKE_BINARY_DIR}/test/magfie/coil_tools_dAphi_dZ_error.dat
        --z-output ${CMAKE_BINARY_DIR}/test/magfie/coil_tools_dAphi_dZ_error.png
        --x-input ${CMAKE_BINARY_DIR}/test/magfie/coil_tools_dAphi_dx_error.dat
        --x-output ${CMAKE_BINARY_DIR}/test/magfie/coil_tools_dAphi_dx_error.png
        --y-input ${CMAKE_BINARY_DIR}/test/magfie/coil_tools_dAphi_dy_error.dat
        --y-output ${CMAKE_BINARY_DIR}/test/magfie/coil_tools_dAphi_dy_error.png
        --dpi 150)
set_tests_properties(plot_vector_potential_deriv_errors PROPERTIES
    DEPENDS test_coil_tools_vector_potential_derivs
    LABELS "magfie;coil_tools;python"
    TIMEOUT 120)

# Register the test
add_test(NAME test_coil_tools_biot_savart COMMAND test_coil_tools_biot_savart.x)

# Set test properties
set_tests_properties(test_coil_tools_biot_savart PROPERTIES 
    LABELS "magfie;biot_savart"
    TIMEOUT 300
)

set(PLOT_FOURIER_SCRIPT ${CMAKE_SOURCE_DIR}/test/magfie/scripts/run_plot_fourier.sh)

add_test(NAME plot_fourier_aug
    COMMAND bash ${PLOT_FOURIER_SCRIPT} ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} aug_lowres)
set_tests_properties(plot_fourier_aug PROPERTIES
    ENVIRONMENT "PYTHON_EXECUTABLE=${Python_EXECUTABLE}"
    LABELS "magfie;coil_tools;python"
    TIMEOUT 600)
