# CMake script to run single-coil Biot-Savart comparison test

set(VACFIELD_EXE "${VACFIELD_EXE}")
set(TEST_DATA_DIR "${TEST_DATA_DIR}")

set(COIL_FILE "${TEST_DATA_DIR}/single_coil.dat")
set(CURRENT_FILE "${TEST_DATA_DIR}/single_coil_currents.txt")
set(GRID_FILE "${TEST_DATA_DIR}/vacfield_single_coil.in")
set(REFERENCE_FILE "${TEST_DATA_DIR}/single_reference.h5")
set(TEST_FILE "${TEST_DATA_DIR}/single_test.nc")
set(PLOT_FILE "${TEST_DATA_DIR}/single_coil_comparison.png")

foreach(path IN LISTS COIL_FILE CURRENT_FILE GRID_FILE)
    if(NOT EXISTS "${path}")
        message(FATAL_ERROR "Required test asset not found: ${path}")
    endif()
endforeach()

file(REMOVE "${REFERENCE_FILE}" "${TEST_FILE}" "${PLOT_FILE}" "${TEST_DATA_DIR}/single_coil_comparison_axis.png")

message(STATUS "Generating reference Fourier data for single coil...")
execute_process(
    COMMAND "${VACFIELD_EXE}" GPEC 1 single_coil.dat
            Fourier vacfield_single_coil.in single_reference.h5
    WORKING_DIRECTORY "${TEST_DATA_DIR}"
    RESULT_VARIABLE result_fourier
    ERROR_VARIABLE error_fourier
)
if(NOT result_fourier EQUAL 0)
    message(FATAL_ERROR "GPEC Fourier failed with exit code ${result_fourier}\n${error_fourier}")
endif()

message(STATUS "Generating vector potential Fourier data for single coil...")
execute_process(
    COMMAND "${VACFIELD_EXE}" GPEC 1 single_coil.dat
            vector_potential vacfield_single_coil.in single_test.nc
    WORKING_DIRECTORY "${TEST_DATA_DIR}"
    RESULT_VARIABLE result_vecpot
    ERROR_VARIABLE error_vecpot
)
if(NOT result_vecpot EQUAL 0)
    message(FATAL_ERROR "vector_potential run failed with exit code ${result_vecpot}\n${error_vecpot}")
endif()

if(NOT EXISTS "${REFERENCE_FILE}" OR NOT EXISTS "${TEST_FILE}")
    message(FATAL_ERROR "Expected output files were not produced")
endif()

if(NOT DEFINED PROJECT_SOURCE_DIR)
    get_filename_component(PROJECT_SOURCE_DIR "${TEST_DATA_DIR}/../../.." ABSOLUTE)
endif()

set(SCRIPT "${PROJECT_SOURCE_DIR}/python/scripts/compare_superposition.py")
if(NOT EXISTS "${SCRIPT}")
    message(FATAL_ERROR "Comparison script not found: ${SCRIPT}")
endif()

message(STATUS "Running comparison for single coil scenario...")
execute_process(
    COMMAND python3 "${SCRIPT}"
            "${REFERENCE_FILE}"
            "${TEST_FILE}"
            "${CURRENT_FILE}"
            "${COIL_FILE}"
            -o "${PLOT_FILE}"
            --ntor 0
            --axis-origin 197.2682697 72.00853957 78.0
            --axis-normal 0.35045589 -0.45058614 0.82106808
            --coil-radius 35.0
            --axis-range 60.0
            --axis-samples 181
    WORKING_DIRECTORY "${TEST_DATA_DIR}"
    RESULT_VARIABLE result_compare
    OUTPUT_VARIABLE output_compare
    ERROR_VARIABLE error_compare
)

if(NOT result_compare EQUAL 0)
    message(FATAL_ERROR "Comparison script failed: ${error_compare}")
endif()

file(WRITE "${TEST_DATA_DIR}/single_coil_summary.txt" "${output_compare}")

message(STATUS "Single coil comparison completed successfully")
message(STATUS "  - ${REFERENCE_FILE}")
message(STATUS "  - ${TEST_FILE}")
message(STATUS "  - ${PLOT_FILE}")
message(STATUS "  - ${TEST_DATA_DIR}/single_coil_comparison_axis.png")
