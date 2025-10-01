set(VACFIELD_EXE "${VACFIELD_EXE}")
set(TEST_DATA_DIR "${TEST_DATA_DIR}")
set(PROJECT_SOURCE_DIR "${PROJECT_SOURCE_DIR}")
set(OUTPUT_DIR "${OUTPUT_DIR}")

if(OUTPUT_DIR STREQUAL "")
    message(FATAL_ERROR "OUTPUT_DIR must be provided")
endif()

set(COIL_FILE "${TEST_DATA_DIR}/axisymmetric_coil.dat")
set(CURRENT_FILE "${TEST_DATA_DIR}/axisymmetric_currents.txt")
set(GRID_FILE "${TEST_DATA_DIR}/vacfield_axisymmetric.in")
set(REFERENCE_FILE "${OUTPUT_DIR}/axisymmetric_reference.h5")
set(TEST_FILE "${OUTPUT_DIR}/axisymmetric_test.nc")
set(SUMMARY_FILE "${OUTPUT_DIR}/axisymmetric_summary.txt")

foreach(path IN LISTS COIL_FILE CURRENT_FILE GRID_FILE)
    if(NOT EXISTS "${path}")
        message(FATAL_ERROR "Required test asset not found: ${path}")
    endif()
endforeach()

file(REMOVE_RECURSE "${OUTPUT_DIR}")
file(MAKE_DIRECTORY "${OUTPUT_DIR}")

message(STATUS "Generating reference Fourier data for axisymmetric coil...")
execute_process(
    COMMAND "${VACFIELD_EXE}" GPEC 1 "${COIL_FILE}"
            Fourier "${GRID_FILE}" "${REFERENCE_FILE}"
    RESULT_VARIABLE result_fourier
    ERROR_VARIABLE error_fourier
)
if(NOT result_fourier EQUAL 0)
    message(FATAL_ERROR "Axisymmetric Fourier run failed with exit code ${result_fourier}\n${error_fourier}")
endif()

message(STATUS "Generating vector potential Fourier data for axisymmetric coil...")
execute_process(
    COMMAND "${VACFIELD_EXE}" GPEC 1 "${COIL_FILE}"
            vector_potential "${GRID_FILE}" "${TEST_FILE}"
    RESULT_VARIABLE result_vecpot
    ERROR_VARIABLE error_vecpot
)
if(NOT result_vecpot EQUAL 0)
    message(FATAL_ERROR "Axisymmetric vector_potential run failed with exit code ${result_vecpot}\n${error_vecpot}")
endif()

if(NOT EXISTS "${REFERENCE_FILE}" OR NOT EXISTS "${TEST_FILE}")
    message(FATAL_ERROR "Expected axisymmetric output files were not produced")
endif()

set(CHECK_SCRIPT "${PROJECT_SOURCE_DIR}/python/tests/check_axisymmetric_loop.py")
if(NOT EXISTS "${CHECK_SCRIPT}")
    message(FATAL_ERROR "Axisymmetric check script not found: ${CHECK_SCRIPT}")
endif()

message(STATUS "Validating axisymmetric coil response...")
execute_process(
    COMMAND python3 "${CHECK_SCRIPT}"
            --reference "${REFERENCE_FILE}"
            --vector "${TEST_FILE}"
            --coil "${COIL_FILE}"
            --currents "${CURRENT_FILE}"
            --summary "${SUMMARY_FILE}"
    WORKING_DIRECTORY "${OUTPUT_DIR}"
    RESULT_VARIABLE result_check
    OUTPUT_VARIABLE output_check
    ERROR_VARIABLE error_check
)

if(NOT result_check EQUAL 0)
    message(FATAL_ERROR "Axisymmetric check failed: ${error_check}")
endif()

message(STATUS "Axisymmetric coil validation completed successfully")
message(STATUS "  - ${REFERENCE_FILE}")
message(STATUS "  - ${TEST_FILE}")
message(STATUS "  - ${SUMMARY_FILE}")
