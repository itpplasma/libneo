# Regenerate a FortSym leaf file into the build tree and require it to match
# the checked-in source. Arguments: GENERATOR, OUTPUT, REFERENCE.
get_filename_component(output_dir "${OUTPUT}" DIRECTORY)
file(MAKE_DIRECTORY "${output_dir}")
execute_process(COMMAND "${GENERATOR}" "${OUTPUT}" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "Generator failed: ${GENERATOR}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${OUTPUT}" "${REFERENCE}"
    RESULT_VARIABLE differ)
if(NOT differ EQUAL 0)
    message(FATAL_ERROR "${REFERENCE} is stale; rebuild the regenerate_ target")
endif()
