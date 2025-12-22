# CheckFortranDependency.cmake
# Provides a function to verify that a Fortran dependency's module files
# are compatible with the current compiler via a smoke test.

# check_fortran_dependency(
#   DEP_NAME        - Name for logging
#   USE_STATEMENT   - Fortran USE statement to test (e.g., "use netcdf")
#   INCLUDE_DIRS    - Include directories for the dependency (module files)
#   LIBRARIES       - Libraries to link (unused but kept for API consistency)
#   RESULT_VAR      - Variable to store result (TRUE if works, FALSE otherwise)
# )
function(check_fortran_dependency DEP_NAME USE_STATEMENT INCLUDE_DIRS LIBRARIES RESULT_VAR)
    # Create test source file
    set(TEST_SOURCE "
program smoke_test
    ${USE_STATEMENT}
    implicit none
    integer :: dummy
    dummy = 0
end program smoke_test
")

    set(TEST_DIR "${CMAKE_BINARY_DIR}/check_${DEP_NAME}")
    file(MAKE_DIRECTORY "${TEST_DIR}")
    set(TEST_FILE "${TEST_DIR}/test.f90")
    file(WRITE "${TEST_FILE}" "${TEST_SOURCE}")

    # Build include flags - also check standard locations
    set(INCLUDE_FLAGS "")

    # Add passed include dirs
    foreach(DIR ${INCLUDE_DIRS})
        if(DIR AND EXISTS "${DIR}")
            list(APPEND INCLUDE_FLAGS "-I${DIR}")
        endif()
    endforeach()

    # Also try /usr/include as fallback for system packages
    if(EXISTS "/usr/include")
        list(APPEND INCLUDE_FLAGS "-I/usr/include")
    endif()

    # Remove duplicates
    list(REMOVE_DUPLICATES INCLUDE_FLAGS)

    # Convert list to space-separated string
    string(REPLACE ";" " " INCLUDE_FLAGS_STR "${INCLUDE_FLAGS}")

    # Try to compile (just syntax check - don't need to link)
    execute_process(
        COMMAND ${CMAKE_Fortran_COMPILER} ${INCLUDE_FLAGS} -fsyntax-only "${TEST_FILE}"
        RESULT_VARIABLE COMPILE_RESULT
        OUTPUT_VARIABLE COMPILE_OUTPUT
        ERROR_VARIABLE COMPILE_ERROR
        TIMEOUT 30
    )

    if(COMPILE_RESULT EQUAL 0)
        message(STATUS "Smoke test for ${DEP_NAME}: PASSED")
        set(${RESULT_VAR} TRUE PARENT_SCOPE)
    else()
        message(STATUS "Smoke test for ${DEP_NAME}: FAILED")
        if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.15")
            message(DEBUG "  Include flags: ${INCLUDE_FLAGS_STR}")
            message(DEBUG "  Error: ${COMPILE_ERROR}")
        endif()
        set(${RESULT_VAR} FALSE PARENT_SCOPE)
    endif()

    # Cleanup
    file(REMOVE_RECURSE "${TEST_DIR}")
endfunction()

# Convenience function to check if we should fetch dependencies
# Returns TRUE if:
#   - FORCE_FETCH_DEPS is ON, or
#   - Compiler is nvfortran/nvhpc (always incompatible with system gfortran libs)
function(should_force_fetch_deps RESULT_VAR)
    if(FORCE_FETCH_DEPS)
        set(${RESULT_VAR} TRUE PARENT_SCOPE)
        return()
    endif()

    # Check for NVIDIA compilers - always need to fetch Fortran-binding deps
    if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR
       CMAKE_Fortran_COMPILER_ID MATCHES "PGI" OR
       CMAKE_Fortran_COMPILER MATCHES "nvfortran")
        set(${RESULT_VAR} TRUE PARENT_SCOPE)
        return()
    endif()

    set(${RESULT_VAR} FALSE PARENT_SCOPE)
endfunction()
