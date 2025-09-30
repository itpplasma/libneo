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

# Generate comparison plot
set(PLOT_OUTPUT "${TEST_DATA_DIR}/comparison.png")

if(DEFINED PYTHON_SCRIPT AND EXISTS "${PYTHON_SCRIPT}")
    message(STATUS "Generating comparison plot...")
    execute_process(
        COMMAND python3 "${PYTHON_SCRIPT}"
                "${TEST_DATA_DIR}/aug_reference.h5"
                "${TEST_DATA_DIR}/aug_test.nc"
                -o "${PLOT_OUTPUT}"
        WORKING_DIRECTORY "${TEST_DATA_DIR}"
        RESULT_VARIABLE result_plot
        OUTPUT_VARIABLE output_plot
        ERROR_VARIABLE error_plot
    )

    if(result_plot EQUAL 0)
        message(STATUS "Comparison plot saved to: ${PLOT_OUTPUT}")
    else()
        message(WARNING "Failed to generate comparison plot: ${error_plot}")
    endif()
else()
    message(WARNING "Comparison script not found: ${PYTHON_SCRIPT}")
endif()

message(STATUS "AUG Biot-Savart comparison test completed successfully")
message(STATUS "Generated files:")
message(STATUS "  - aug_reference.h5 (GPEC Fourier)")
message(STATUS "  - aug_test.nc (coil_tools vector_potential)")
message(STATUS "  - comparison.png (visual comparison)")
