# CMake script to run AUG Biot-Savart comparison test
# This test compares GPEC Fourier vs coil_tools vector_potential implementations

set(VACFIELD_EXE "${VACFIELD_EXE}")
set(TEST_DATA_DIR "${TEST_DATA_DIR}")

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

# Clean up old output files
file(REMOVE "${TEST_DATA_DIR}/aug_reference.h5")
file(REMOVE "${TEST_DATA_DIR}/aug_test.nc")

message(STATUS "Running GPEC Fourier Biot-Savart (aug_reference.h5)...")
execute_process(
    COMMAND "${VACFIELD_EXE}" GPEC 2 aug_bu.dat aug_bl.dat
            Fourier vacfield_AUG_lowres.in aug_reference.h5
    WORKING_DIRECTORY "${TEST_DATA_DIR}"
    RESULT_VARIABLE result_fourier
    OUTPUT_VARIABLE output_fourier
    ERROR_VARIABLE error_fourier
)

if(NOT result_fourier EQUAL 0)
    message(FATAL_ERROR "GPEC Fourier failed with exit code ${result_fourier}\n${error_fourier}")
endif()

message(STATUS "Running coil_tools vector_potential Biot-Savart (aug_test.nc)...")
execute_process(
    COMMAND "${VACFIELD_EXE}" GPEC 2 aug_bu.dat aug_bl.dat
            vector_potential vacfield_AUG_lowres.in aug_test.nc
    WORKING_DIRECTORY "${TEST_DATA_DIR}"
    RESULT_VARIABLE result_vecpot
    OUTPUT_VARIABLE output_vecpot
    ERROR_VARIABLE error_vecpot
)

if(NOT result_vecpot EQUAL 0)
    message(FATAL_ERROR "vector_potential failed with exit code ${result_vecpot}\n${error_vecpot}")
endif()

# Check that output files were created
if(NOT EXISTS "${TEST_DATA_DIR}/aug_reference.h5")
    message(FATAL_ERROR "Output file aug_reference.h5 was not created")
endif()

if(NOT EXISTS "${TEST_DATA_DIR}/aug_test.nc")
    message(FATAL_ERROR "Output file aug_test.nc was not created")
endif()

message(STATUS "Generated field files successfully")

# Generate superposition comparison plot and check tolerance
if(NOT DEFINED PROJECT_SOURCE_DIR)
    get_filename_component(PROJECT_SOURCE_DIR "${TEST_DATA_DIR}/../../.." ABSOLUTE)
endif()

set(SUPERPOSITION_SCRIPT "${PROJECT_SOURCE_DIR}/python/scripts/compare_superposition.py")
set(SUPERPOSITION_PLOT "${TEST_DATA_DIR}/superposition_comparison.png")

if(EXISTS "${SUPERPOSITION_SCRIPT}")
    message(STATUS "Generating three-way comparison and checking tolerance...")
    execute_process(
        COMMAND python3 "${SUPERPOSITION_SCRIPT}"
                "${TEST_DATA_DIR}/aug_reference.h5"
                "${TEST_DATA_DIR}/aug_test.nc"
                "${TEST_DATA_DIR}/aug_currents.txt"
                "${TEST_DATA_DIR}/aug_bu.dat"
                "${TEST_DATA_DIR}/aug_bl.dat"
                -o "${SUPERPOSITION_PLOT}"
        WORKING_DIRECTORY "${TEST_DATA_DIR}"
        RESULT_VARIABLE result_superpos
        OUTPUT_VARIABLE output_superpos
        ERROR_VARIABLE error_superpos
    )

    if(NOT result_superpos EQUAL 0)
        message(FATAL_ERROR "Comparison failed: ${error_superpos}")
    endif()

    # Parse output to extract median and mean relative errors
    string(REGEX MATCH "Median relative error: ([0-9.]+)%" median_match "${output_superpos}")
    string(REGEX MATCH "Mean relative error: +([0-9.]+)%" mean_match "${output_superpos}")

    if(median_match)
        string(REGEX REPLACE "Median relative error: ([0-9.]+)%" "\\1" median_error "${median_match}")
        message(STATUS "Median relative error: ${median_error}%")

        # Check tolerance (accept up to 50% median error)
        if(median_error GREATER 50.0)
            message(FATAL_ERROR "Median relative error (${median_error}%) exceeds tolerance (50%)")
        endif()
    endif()

    if(mean_match)
        string(REGEX REPLACE "Mean relative error: +([0-9.]+)%" "\\1" mean_error "${mean_match}")
        message(STATUS "Mean relative error: ${mean_error}%")

        # Check tolerance (accept up to 70% mean error)
        if(mean_error GREATER 70.0)
            message(FATAL_ERROR "Mean relative error (${mean_error}%) exceeds tolerance (70%)")
        endif()
    endif()

    message(STATUS "Comparison plot saved to: ${SUPERPOSITION_PLOT}")
else()
    message(FATAL_ERROR "Comparison script not found: ${SUPERPOSITION_SCRIPT}")
endif()

message(STATUS "AUG Biot-Savart comparison test completed successfully")
message(STATUS "Generated files:")
message(STATUS "  - aug_reference.h5 (GPEC Fourier)")
message(STATUS "  - aug_test.nc (coil_tools vector_potential)")
message(STATUS "  - superposition_comparison.png (comparison plot)")
