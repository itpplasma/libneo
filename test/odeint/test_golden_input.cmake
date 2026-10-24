set(input "${TEST_DIRECTORY}/reference.f90")
set(output "${TEST_DIRECTORY}/generated/reference.f90")
file(MAKE_DIRECTORY "${TEST_DIRECTORY}")
file(WRITE "${input}" "module public_reference\nend module public_reference\n")
file(SHA256 "${input}" expected_sha256)

execute_process(
    COMMAND "${CMAKE_COMMAND}"
        "-DINPUT=${input}" "-DOUTPUT=${output}"
        "-DEXPECTED_SHA256=${expected_sha256}" -P "${PREPARE_SCRIPT}"
    RESULT_VARIABLE result
)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Valid immutable golden input was rejected")
endif()
file(SHA256 "${output}" actual_sha256)
if(NOT actual_sha256 STREQUAL expected_sha256)
    message(FATAL_ERROR "Immutable golden bytes changed while copying")
endif()

file(WRITE "${input}" "module changed_reference\nend module changed_reference\n")
execute_process(
    COMMAND "${CMAKE_COMMAND}"
        "-DINPUT=${input}" "-DOUTPUT=${output}"
        "-DEXPECTED_SHA256=${expected_sha256}" -P "${PREPARE_SCRIPT}"
    RESULT_VARIABLE result OUTPUT_QUIET ERROR_QUIET
)
if(result EQUAL 0)
    message(FATAL_ERROR "Changed immutable golden input was accepted")
endif()
file(SHA256 "${output}" actual_sha256)
if(NOT actual_sha256 STREQUAL expected_sha256)
    message(FATAL_ERROR "Rejected input replaced the verified golden output")
endif()
