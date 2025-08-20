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

# Golden record test - dynamically fetches and builds reference version
add_executable(test_golden_record_odeint.x 
    test_golden_record_odeint.f90
    golden_odeint_dynamic.f90
)
target_link_libraries(test_golden_record_odeint.x ${COMMON_LIBS} odeint_test_common)

# Add a pre-build step to fetch golden record
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/golden_odeint_dynamic.f90
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/fetch_golden_record.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Fetching golden record from latest tagged version"
)

# Ensure the golden file is built before the test
add_custom_target(golden_fetch_odeint
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/golden_odeint_dynamic.f90
)
add_dependencies(test_golden_record_odeint.x golden_fetch_odeint)

# Register tests
function (add_fortran_test name)
  add_test(NAME ${name} COMMAND ${name}.x)
  set_tests_properties(${name} PROPERTIES  FAIL_REGULAR_EXPRESSION "STOP")
endfunction()

add_fortran_test(test_odeint_allroutines_context)
add_fortran_test(test_odeint_thread_safety) 
add_fortran_test(test_golden_record_odeint)