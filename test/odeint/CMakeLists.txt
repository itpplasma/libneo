set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Common test utilities library
add_library(odeint_test_common odeint_test_common.f90)
target_link_libraries(odeint_test_common ${COMMON_LIBS})

# Serial context test
add_executable(test_odeint_allroutines_context.x test_odeint_allroutines_context.f90)
target_link_libraries(test_odeint_allroutines_context.x ${COMMON_LIBS} odeint_test_common)

# Thread safety test
add_executable(test_odeint_thread_safety.x test_odeint_thread_safety.f90)
target_link_libraries(test_odeint_thread_safety.x ${COMMON_LIBS} odeint_test_common)
target_compile_options(test_odeint_thread_safety.x PRIVATE ${OpenMP_Fortran_FLAGS})
target_link_libraries(test_odeint_thread_safety.x ${OpenMP_Fortran_LIBRARIES})

if(LIBNEO_ENABLE_GOLDEN_TESTS)
    set(LIBNEO_ODEINT_GOLDEN_SOURCE_DIR "" CACHE PATH
        "Directory containing an immutable golden_odeint_dynamic.f90")
    set(LIBNEO_ODEINT_GOLDEN_SHA256 "" CACHE STRING
        "Expected SHA256 of the immutable ODEINT golden source")
    set(_golden_output "${CMAKE_CURRENT_BINARY_DIR}/golden_odeint_dynamic.f90")
    if(LIBNEO_ODEINT_GOLDEN_SOURCE_DIR)
        set(_golden_input
            "${LIBNEO_ODEINT_GOLDEN_SOURCE_DIR}/golden_odeint_dynamic.f90")
        file(SHA256 "${_golden_input}" _golden_sha256)
        if(NOT _golden_sha256 STREQUAL LIBNEO_ODEINT_GOLDEN_SHA256)
            message(FATAL_ERROR "Immutable ODEINT golden source SHA256 mismatch")
        endif()
        add_custom_command(
            OUTPUT "${_golden_output}"
            COMMAND ${CMAKE_COMMAND}
                "-DINPUT=${_golden_input}" "-DOUTPUT=${_golden_output}"
                "-DEXPECTED_SHA256=${LIBNEO_ODEINT_GOLDEN_SHA256}"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/prepare_golden_input.cmake"
            DEPENDS "${_golden_input}"
                "${CMAKE_CURRENT_SOURCE_DIR}/prepare_golden_input.cmake"
            COMMENT "Verifying immutable ODEINT golden source"
            VERBATIM
        )
    else()
        if(LIBNEO_ODEINT_GOLDEN_SHA256)
            message(FATAL_ERROR "ODEINT golden SHA256 requires a source directory")
        endif()
        add_custom_command(
            OUTPUT "${_golden_output}"
            COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/fetch_golden_record.sh"
                "${_golden_output}"
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            COMMENT "Fetching golden record from latest tagged version"
            VERBATIM
        )
    endif()

    add_executable(test_golden_record_odeint.x EXCLUDE_FROM_ALL
        test_golden_record_odeint.f90
        "${_golden_output}"
    )
    target_link_libraries(test_golden_record_odeint.x ${COMMON_LIBS} odeint_test_common)
    if(LIBNEO_ODEINT_GOLDEN_SOURCE_DIR)
        # Name the provider input in the native File API as well as the recipe.
        set_source_files_properties("${_golden_input}" PROPERTIES HEADER_FILE_ONLY TRUE)
        target_sources(test_golden_record_odeint.x PRIVATE "${_golden_input}")
    endif()

    # Ensure the golden file is built before the test
    add_custom_target(golden_fetch_odeint
        DEPENDS "${_golden_output}"
    )
    add_dependencies(test_golden_record_odeint.x golden_fetch_odeint)
else()
    message(STATUS
        "Skipping ODEINT golden record regression test; set LIBNEO_ENABLE_GOLDEN_TESTS=ON to enable.")
endif()

# Register tests
add_test(NAME test_odeint_allroutines_context COMMAND test_odeint_allroutines_context.x)
add_test(NAME test_odeint_thread_safety COMMAND test_odeint_thread_safety.x)

if(LIBNEO_ENABLE_GOLDEN_TESTS)
    add_test(NAME test_odeint_golden_input_contract
        COMMAND ${CMAKE_COMMAND}
            "-DTEST_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}/golden-input-contract"
            "-DPREPARE_SCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/prepare_golden_input.cmake"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/test_golden_input.cmake"
    )
    add_test(NAME build_test_golden_record_odeint
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
        --target test_golden_record_odeint.x
    )
    set_tests_properties(build_test_golden_record_odeint PROPERTIES
        FIXTURES_SETUP odeint_golden_fixture
        LABELS build-helper
    )

    add_test(NAME test_golden_record_odeint COMMAND test_golden_record_odeint.x)
    set_tests_properties(test_golden_record_odeint PROPERTIES
        FIXTURES_REQUIRED odeint_golden_fixture
    )
endif()

# Event detection test
add_executable(test_odeint_events.x test_odeint_events.f90)
target_link_libraries(test_odeint_events.x ${COMMON_LIBS})
add_test(NAME test_odeint_events COMMAND test_odeint_events.x)

add_executable(test_odeint_scale_floor.x test_odeint_scale_floor.f90)
target_link_libraries(test_odeint_scale_floor.x ${COMMON_LIBS})
add_test(NAME test_odeint_scale_floor COMMAND test_odeint_scale_floor.x)
