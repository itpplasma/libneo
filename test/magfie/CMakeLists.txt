set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/magfie)

if(NOT Python_FOUND)
    message(FATAL_ERROR
        "Python support is required for magfie tests; ensure libneo top-level "
        "find_package(Python ...) succeeds before configuring this directory."
    )
endif()

set(_MAGFIE_TEST_PYTHON_EXECUTABLE "${Python_EXECUTABLE}")

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

# Register the test
add_test(NAME test_coil_tools_biot_savart COMMAND test_coil_tools_biot_savart.x)

# Set test properties
set_tests_properties(test_coil_tools_biot_savart PROPERTIES
    LABELS "magfie;biot_savart"
    TIMEOUT 300
)

add_executable(test_gpec_loader_comparison.x test_gpec_loader_comparison.f90)
target_link_libraries(test_gpec_loader_comparison.x PRIVATE
    neo
    magfie
    ${FFTW3_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
)

add_test(NAME test_gpec_loader_comparison
    COMMAND test_gpec_loader_comparison.x ${CMAKE_SOURCE_DIR}/test/magfie/test_data
)

set_tests_properties(test_gpec_loader_comparison PROPERTIES
    LABELS "magfie;gpec;loader"
    TIMEOUT 120
)

# Integration test: Compare GPEC Fourier vs coil_tools vector_potential on AUG data
add_test(NAME test_aug_biot_savart_comparison
    COMMAND ${CMAKE_COMMAND}
        -DVACFIELD_EXE=$<TARGET_FILE:vacfield.x>
        -DTEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/test_data
        -DPROJECT_SOURCE_DIR=${CMAKE_SOURCE_DIR}
        -DPYTHON_EXECUTABLE=${_MAGFIE_TEST_PYTHON_EXECUTABLE}
        -DOUTPUT_DIR:PATH=${CMAKE_BINARY_DIR}/test/magfie
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
    fortplot::fortplot
    ${FFTW3_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
)

add_test(NAME benchmark_biot_savart
    COMMAND ${CMAKE_COMMAND}
        -DBENCHMARK_EXE=$<TARGET_FILE:benchmark_biot_savart.x>
        -DTEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/test_data
        -DOUTPUT_DIR=${CMAKE_BINARY_DIR}/test/magfie/benchmark
        -P ${CMAKE_CURRENT_SOURCE_DIR}/run_benchmark_case.cmake
)

set_tests_properties(benchmark_biot_savart PROPERTIES
    LABELS "magfie;biot_savart;benchmark"
    DEPENDS benchmark_biot_savart.x
)

# Integration test: single-coil scenario with analytical verification
add_test(NAME test_single_coil_biot_savart
    COMMAND ${CMAKE_COMMAND}
        -DVACFIELD_EXE=$<TARGET_FILE:vacfield.x>
        -DTEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/test_data
        -DPROJECT_SOURCE_DIR=${CMAKE_SOURCE_DIR}
        -DPYTHON_EXECUTABLE=${_MAGFIE_TEST_PYTHON_EXECUTABLE}
        -DOUTPUT_DIR:PATH=${CMAKE_BINARY_DIR}/test/magfie/single_coil
        -P ${CMAKE_CURRENT_SOURCE_DIR}/run_single_coil_case.cmake
)

set_tests_properties(test_single_coil_biot_savart PROPERTIES
    LABELS "magfie;biot_savart;integration;single_coil"
    TIMEOUT 600
    DEPENDS vacfield.x
)

add_test(NAME test_axisymmetric_biot_savart
    COMMAND ${CMAKE_COMMAND}
        -DVACFIELD_EXE=$<TARGET_FILE:vacfield.x>
        -DTEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/test_data
        -DPROJECT_SOURCE_DIR=${CMAKE_SOURCE_DIR}
        -DPYTHON_EXECUTABLE=${_MAGFIE_TEST_PYTHON_EXECUTABLE}
        -DOUTPUT_DIR:PATH=${CMAKE_BINARY_DIR}/test/magfie/axisymmetric_coil
        -P ${CMAKE_CURRENT_SOURCE_DIR}/run_axisymmetric_case.cmake
)

set_tests_properties(test_axisymmetric_biot_savart PROPERTIES
    LABELS "magfie;biot_savart;integration;axisymmetric"
    TIMEOUT 300
    DEPENDS vacfield.x
)


add_executable(test_coil_tools_segment_formulas.x
    test_coil_tools_segment_formulas.f90)
target_link_libraries(test_coil_tools_segment_formulas.x PRIVATE
    neo
    magfie
    util_for_test
    ${FFTW3_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
)
target_compile_definitions(test_coil_tools_segment_formulas.x PRIVATE
    COIL_TOOLS_PLOT_DIR="${CMAKE_BINARY_DIR}/test/magfie")

add_test(NAME test_coil_tools_segment_formulas
    COMMAND test_coil_tools_segment_formulas.x)
set_tests_properties(test_coil_tools_segment_formulas PROPERTIES
    LABELS "magfie;coil_tools;formulas"
    TIMEOUT 60
)

add_executable(test_coil_tools_aphi_gradients.x
    test_coil_tools_aphi_gradients.f90)
target_link_libraries(test_coil_tools_aphi_gradients.x PRIVATE
    neo
    magfie
    util_for_test
    ${FFTW3_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
)
target_compile_definitions(test_coil_tools_aphi_gradients.x PRIVATE
    COIL_TOOLS_PLOT_DIR="${CMAKE_BINARY_DIR}/test/magfie"
    COIL_TOOLS_SCRIPT_DIR="${CMAKE_SOURCE_DIR}/test/magfie")

add_test(NAME test_coil_tools_aphi_gradients
    COMMAND test_coil_tools_aphi_gradients.x)
set_tests_properties(test_coil_tools_aphi_gradients PROPERTIES
    LABELS "magfie;coil_tools;derivatives"
    TIMEOUT 60
)
