set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/magfie)

include(FetchContent)

set(FORTPLOT_OVERRIDE_REPO $ENV{FORTPLOT_OVERRIDE_REPO})
if (FORTPLOT_OVERRIDE_REPO)
    message(STATUS "Using overridden fortplot repository: ${FORTPLOT_OVERRIDE_REPO}")
    set(_fortplot_repo ${FORTPLOT_OVERRIDE_REPO})
    set(_fortplot_tag HEAD)
else()
    set(_fortplot_repo https://github.com/lazy-fortran/fortplot.git)
    set(_fortplot_tag main)
endif()

FetchContent_Declare(fortplot
    GIT_REPOSITORY ${_fortplot_repo}
    GIT_TAG ${_fortplot_tag}
    GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(fortplot)

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
    fortplot::fortplot
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

# Register the test
add_test(NAME test_coil_tools_biot_savart COMMAND test_coil_tools_biot_savart.x)

# Set test properties
set_tests_properties(test_coil_tools_biot_savart PROPERTIES
    LABELS "magfie;biot_savart"
    TIMEOUT 300
)

# Test for ntor=0 field reconstruction with tilted coil
add_executable(test_ntor0_tilted_coil.x test_ntor0_tilted_coil.f90)
target_link_libraries(test_ntor0_tilted_coil.x PRIVATE
    neo
    magfie
    util_for_test
    ${FFTW3_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
)

add_executable(compute_axis_biot_savart.x compute_axis_biot_savart.f90)
target_link_libraries(compute_axis_biot_savart.x PRIVATE
    neo
    magfie
)

add_test(
    NAME generate_tilted_coil_geometry
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generate_tilted_coil.py
            ${CMAKE_CURRENT_BINARY_DIR}/tilted_coil.dat
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(
    NAME compute_tilted_coil_fourier_modes
    COMMAND test_ntor0_tilted_coil.x
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(compute_tilted_coil_fourier_modes PROPERTIES
    DEPENDS generate_tilted_coil_geometry
    LABELS "magfie;ntor0"
    TIMEOUT 300
)

add_test(
    NAME validate_ntor0_field_reconstruction
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/validate_ntor0_fix.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(validate_ntor0_field_reconstruction PROPERTIES
    DEPENDS compute_tilted_coil_fourier_modes
    LABELS "magfie;ntor0"
    TIMEOUT 300
    WILL_FAIL TRUE
)
