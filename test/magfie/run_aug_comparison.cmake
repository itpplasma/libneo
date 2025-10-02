# CMake script to run AUG Biot-Savart comparison test
# This test compares GPEC Fourier vs coil_tools vector_potential implementations

set(VACFIELD_EXE "${VACFIELD_EXE}")
set(TEST_DATA_DIR "${TEST_DATA_DIR}")
set(OUTPUT_DIR "${OUTPUT_DIR}")
set(AUG_OUTPUT_DIR "${OUTPUT_DIR}/aug")

# Check that required files exist
set(REQUIRED_FILES
    "${TEST_DATA_DIR}/aug_bu.dat"
    "${TEST_DATA_DIR}/aug_bl.dat"
    "${TEST_DATA_DIR}/vacfield_AUG_lowres.in"
)

foreach(file ${REQUIRED_FILES})
    if(NOT EXISTS "${file}")
        message(FATAL_ERROR "Required test file not found: ${file}")
    endif()
endforeach()

# Create output directory hierarchy
file(MAKE_DIRECTORY "${AUG_OUTPUT_DIR}")

# Clean up old output files
file(REMOVE "${AUG_OUTPUT_DIR}/aug_reference.h5")
file(REMOVE "${AUG_OUTPUT_DIR}/aug_test.nc")

message(STATUS "Running GPEC Fourier Biot-Savart (aug_reference.h5)...")
execute_process(
    COMMAND "${VACFIELD_EXE}" GPEC 2 ${TEST_DATA_DIR}/aug_bu.dat ${TEST_DATA_DIR}/aug_bl.dat
            Fourier ${TEST_DATA_DIR}/vacfield_AUG_lowres.in aug_reference.h5
    WORKING_DIRECTORY "${AUG_OUTPUT_DIR}"
    RESULT_VARIABLE result_fourier
    OUTPUT_VARIABLE output_fourier
    ERROR_VARIABLE error_fourier
)

if(NOT result_fourier EQUAL 0)
    message(FATAL_ERROR "GPEC Fourier failed with exit code ${result_fourier}\n${error_fourier}")
endif()

message(STATUS "Running coil_tools vector_potential Biot-Savart (aug_test.nc)...")
execute_process(
    COMMAND "${VACFIELD_EXE}" GPEC 2 ${TEST_DATA_DIR}/aug_bu.dat ${TEST_DATA_DIR}/aug_bl.dat
            vector_potential ${TEST_DATA_DIR}/vacfield_AUG_lowres.in aug_test.nc
    WORKING_DIRECTORY "${AUG_OUTPUT_DIR}"
    RESULT_VARIABLE result_vecpot
    OUTPUT_VARIABLE output_vecpot
    ERROR_VARIABLE error_vecpot
)

if(NOT result_vecpot EQUAL 0)
    message(FATAL_ERROR "vector_potential failed with exit code ${result_vecpot}\n${error_vecpot}")
endif()

# Check that output files were created
if(NOT EXISTS "${AUG_OUTPUT_DIR}/aug_reference.h5")
    message(FATAL_ERROR "Output file aug_reference.h5 was not created")
endif()

if(NOT EXISTS "${AUG_OUTPUT_DIR}/aug_test.nc")
    message(FATAL_ERROR "Output file aug_test.nc was not created")
endif()

message(STATUS "Generated field files successfully")

# Generate superposition comparison plot and check tolerance
if(NOT DEFINED PROJECT_SOURCE_DIR)
    get_filename_component(PROJECT_SOURCE_DIR "${TEST_DATA_DIR}/../../.." ABSOLUTE)
endif()

set(COMPARISON_SCRIPT "${PROJECT_SOURCE_DIR}/python/scripts/plot_biotsavart_fourier.py")
set(PER_COIL_PLOT "${AUG_OUTPUT_DIR}/biotsavart_fourier.png")
set(SUM_PLOT "${AUG_OUTPUT_DIR}/biotsavart_fourier_sum.png")

set(MEDIAN_ERROR_EXCEEDED FALSE)
set(MEAN_ERROR_EXCEEDED FALSE)
set(MEDIAN_ERROR_VALUE 0.0)
set(MEAN_ERROR_VALUE 0.0)

if(NOT DEFINED PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE python3)
endif()

if(EXISTS "${COMPARISON_SCRIPT}")
    message(STATUS "Running unified Biot-Savart comparison...")
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env MPLBACKEND=Agg
                "${PYTHON_EXECUTABLE}" "${COMPARISON_SCRIPT}"
                aug_reference.h5
                aug_test.nc
                --currents "${TEST_DATA_DIR}/aug_currents.txt"
                --coil-files "${TEST_DATA_DIR}/aug_bu.dat" "${TEST_DATA_DIR}/aug_bl.dat"
                --ntor 2
                --per-coil-output biotsavart_fourier.png
                --sum-output biotsavart_fourier_sum.png
                --deriv-diff-output biotsavart_fourier_derivdiff.png
        WORKING_DIRECTORY "${AUG_OUTPUT_DIR}"
        RESULT_VARIABLE comparison_result
        OUTPUT_VARIABLE comparison_output
        ERROR_VARIABLE comparison_error
    )

    if(NOT comparison_result EQUAL 0)
        message(FATAL_ERROR "Comparison script failed: ${comparison_error}")
    endif()

    string(REGEX MATCH "Median relative error: ([0-9.]+)%" median_match "${comparison_output}")
    string(REGEX MATCH "Mean relative error: +([0-9.]+)%" mean_match "${comparison_output}")

    if(median_match)
        string(REGEX REPLACE "Median relative error: ([0-9.]+)%" "\\1" MEDIAN_ERROR_VALUE "${median_match}")
        message(STATUS "Median relative error: ${MEDIAN_ERROR_VALUE}%")
        if(MEDIAN_ERROR_VALUE GREATER 50.0)
            set(MEDIAN_ERROR_EXCEEDED TRUE)
        endif()
    endif()

    if(mean_match)
        string(REGEX REPLACE "Mean relative error: +([0-9.]+)%" "\\1" MEAN_ERROR_VALUE "${mean_match}")
        message(STATUS "Mean relative error: ${MEAN_ERROR_VALUE}%")
        if(MEAN_ERROR_VALUE GREATER 70.0)
            set(MEAN_ERROR_EXCEEDED TRUE)
        endif()
    endif()

    message(STATUS "Per-coil plot saved to: ${PER_COIL_PLOT}")
    message(STATUS "Summed plot saved to: ${SUM_PLOT}")
else()
    message(FATAL_ERROR "Comparison script not found: ${COMPARISON_SCRIPT}")
endif()

message(STATUS "AUG Biot-Savart comparison test completed successfully")
message(STATUS "Generated files under ${AUG_OUTPUT_DIR}:")
message(STATUS "  - aug_reference.h5 (GPEC Fourier)")
message(STATUS "  - aug_test.nc (coil_tools vector_potential)")
message(STATUS "  - biotsavart_fourier.png (per-coil diagnostics)")
message(STATUS "  - biotsavart_fourier_sum.png (summed-field diagnostics)")
message(STATUS "  - biotsavart_fourier_derivdiff.png (gauge derivative diagnostics)")

if(MEDIAN_ERROR_EXCEEDED)
    message(FATAL_ERROR "Median relative error (${MEDIAN_ERROR_VALUE}%) exceeds tolerance (50%)")
endif()

if(MEAN_ERROR_EXCEEDED)
    message(FATAL_ERROR "Mean relative error (${MEAN_ERROR_VALUE}%) exceeds tolerance (70%)")
endif()
