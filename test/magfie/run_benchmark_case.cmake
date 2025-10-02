set(BENCHMARK_EXE "${BENCHMARK_EXE}")
set(TEST_DATA_DIR "${TEST_DATA_DIR}")
set(OUTPUT_DIR "${OUTPUT_DIR}")

if(OUTPUT_DIR STREQUAL "")
    message(FATAL_ERROR "OUTPUT_DIR must be provided")
endif()

if(NOT EXISTS "${BENCHMARK_EXE}")
    message(FATAL_ERROR "Benchmark executable not found: ${BENCHMARK_EXE}")
endif()

if(NOT IS_DIRECTORY "${TEST_DATA_DIR}")
    message(FATAL_ERROR "Test data directory not found: ${TEST_DATA_DIR}")
endif()

file(MAKE_DIRECTORY "${OUTPUT_DIR}")

set(SUMMARY_FILE "${OUTPUT_DIR}/benchmark_summary.txt")
set(ERROR_PLOT "${OUTPUT_DIR}/benchmark_biot_savart_errors.png")

execute_process(
    COMMAND "${BENCHMARK_EXE}" "${TEST_DATA_DIR}" "${OUTPUT_DIR}"
    RESULT_VARIABLE result_bench
    OUTPUT_VARIABLE output_bench
    ERROR_VARIABLE error_bench
    WORKING_DIRECTORY "${OUTPUT_DIR}"
)

file(WRITE "${SUMMARY_FILE}" "${output_bench}")

if(NOT result_bench EQUAL 0)
    message(FATAL_ERROR "benchmark_biot_savart failed: ${error_bench}")
endif()

if(NOT EXISTS "${ERROR_PLOT}")
    message(FATAL_ERROR "Benchmark did not produce ${ERROR_PLOT}")
endif()

message(STATUS "Benchmark outputs located in ${OUTPUT_DIR}")
