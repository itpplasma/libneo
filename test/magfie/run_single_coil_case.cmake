# CMake script to run single-coil Biot-Savart comparison test

set(VACFIELD_EXE "${VACFIELD_EXE}")
set(TEST_DATA_DIR "${TEST_DATA_DIR}")
set(OUTPUT_DIR "${OUTPUT_DIR}")

if(OUTPUT_DIR STREQUAL "")
    message(FATAL_ERROR "OUTPUT_DIR must be provided")
endif()

set(COIL_FILE "${TEST_DATA_DIR}/single_coil.dat")
set(CURRENT_FILE "${TEST_DATA_DIR}/single_coil_currents.txt")
set(GRID_FILE "${TEST_DATA_DIR}/vacfield_single_coil.in")
set(REFERENCE_FILE "${OUTPUT_DIR}/single_reference.h5")
set(TEST_FILE "${OUTPUT_DIR}/single_test.nc")
set(PLOT_FILE "${OUTPUT_DIR}/single_coil_per_coil.png")
set(SUM_PLOT "${OUTPUT_DIR}/single_coil_sum.png")
set(AXIS_PLOT "${OUTPUT_DIR}/single_coil_comparison_axis.png")
set(SUMMARY_FILE "${OUTPUT_DIR}/single_coil_summary.txt")

foreach(path IN LISTS COIL_FILE CURRENT_FILE GRID_FILE)
    if(NOT EXISTS "${path}")
        message(FATAL_ERROR "Required test asset not found: ${path}")
    endif()
endforeach()

file(REMOVE_RECURSE "${OUTPUT_DIR}")
file(MAKE_DIRECTORY "${OUTPUT_DIR}")

file(REMOVE
    "${REFERENCE_FILE}"
    "${TEST_FILE}"
    "${PLOT_FILE}"
    "${SUM_PLOT}"
    "${AXIS_PLOT}"
    "${SUMMARY_FILE}"
)

message(STATUS "Generating reference Fourier data for single coil...")
execute_process(
    COMMAND "${VACFIELD_EXE}" GPEC 1 "${COIL_FILE}"
            Fourier "${GRID_FILE}" "${REFERENCE_FILE}"
    RESULT_VARIABLE result_fourier
    ERROR_VARIABLE error_fourier
)
if(NOT result_fourier EQUAL 0)
    message(FATAL_ERROR "GPEC Fourier failed with exit code ${result_fourier}\n${error_fourier}")
endif()

message(STATUS "Generating vector potential Fourier data for single coil...")
execute_process(
    COMMAND "${VACFIELD_EXE}" GPEC 1 "${COIL_FILE}"
            vector_potential "${GRID_FILE}" "${TEST_FILE}"
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

set(COMPARISON_SCRIPT "${PROJECT_SOURCE_DIR}/python/scripts/plot_biotsavart_fourier.py")
if(NOT EXISTS "${COMPARISON_SCRIPT}")
    message(FATAL_ERROR "Comparison script not found: ${COMPARISON_SCRIPT}")
endif()

message(STATUS "Running comparison for single coil scenario...")
execute_process(
    COMMAND python3 "${COMPARISON_SCRIPT}"
            "${REFERENCE_FILE}"
            "${TEST_FILE}"
            --currents "${CURRENT_FILE}"
            --coil-files "${COIL_FILE}"
            --ntor 0
            --per-coil-output "${PLOT_FILE}"
            --sum-output "${SUM_PLOT}"
            --axis-output "${AXIS_PLOT}"
            --axis-origin 197.2682697 72.00853957 78.0
            --axis-normal 0.35045589 -0.45058614 0.82106808
            --coil-radius 35.0
            --axis-range 60.0
            --axis-samples 181
    WORKING_DIRECTORY "${OUTPUT_DIR}"
    RESULT_VARIABLE result_compare
    OUTPUT_VARIABLE output_compare
    ERROR_VARIABLE error_compare
)

if(NOT result_compare EQUAL 0)
    message(FATAL_ERROR "Comparison script failed: ${error_compare}")
endif()

file(WRITE "${SUMMARY_FILE}" "${output_compare}")

message(STATUS "Single coil comparison completed successfully")
message(STATUS "  - ${REFERENCE_FILE}")
message(STATUS "  - ${TEST_FILE}")
message(STATUS "  - ${PLOT_FILE}")
message(STATUS "  - ${SUM_PLOT}")
message(STATUS "  - ${AXIS_PLOT}")
