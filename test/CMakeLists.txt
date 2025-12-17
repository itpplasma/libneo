set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set (MAIN_LIB ${MAIN_LIB} neo)
set (HDF5_TOOLS_LIB ${HDF5_TOOLS_LIB} hdf5_tools)
set (MAGFIE_LIB ${MAGFIE_LIB} magfie)

set(PROJLIBS /proj/plasma/Libs/ CACHE STRING "Common library path")
include_directories(${MPI_Fortran_INCLUDE_PATH})

include_directories(${PROJECT_SOURCE_DIR}/src/hdf5_tools)
include_directories(${PROJECT_SOURCE_DIR}/extra/MyMPILib)
include_directories(${MPI_Fortran_INCLUDE_PATH})

# Assumes previous include of hdf5_tools/CMakeLists.txt for find_package
include_directories(${HDF5_INCLUDE_DIRS})

set(LIBNEO_TESTING_ENABLED OFF)
if(DEFINED ENV{LIBNEO_TESTING})
  string(STRIP "$ENV{LIBNEO_TESTING}" _libneo_testing_raw)
  string(TOLOWER "${_libneo_testing_raw}" _libneo_testing_value)
  if(NOT _libneo_testing_value STREQUAL ""
      AND NOT _libneo_testing_value STREQUAL "0"
      AND NOT _libneo_testing_value STREQUAL "false"
      AND NOT _libneo_testing_value STREQUAL "no"
      AND NOT _libneo_testing_value STREQUAL "off")
    set(LIBNEO_TESTING_ENABLED ON)
  endif()
endif()

find_package(Python COMPONENTS Interpreter)
if (NOT Python_FOUND)
  message(FATAL_ERROR "Python interpreter required for GEQDSK scripts")
endif()

execute_process(COMMAND nf-config --flibs
                OUTPUT_VARIABLE NETCDF_LIBRARIES
                OUTPUT_STRIP_TRAILING_WHITESPACE)

set (COMMON_LIBS ${MAIN_LIB} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES})

add_executable(test_transport.x transport/test_transport.f90)
target_link_libraries(test_transport.x
  neo
  transport
  collision_freqs
  util_for_test
)

add_executable(test_collision_freqs.x collisions/test_collision_freqs.f90)
target_link_libraries(test_collision_freqs.x
  neo
  collision_freqs
  species
  util_for_test
)

add_executable(test_arnoldi.x source/test_arnoldi.f90)
target_link_libraries(test_arnoldi.x
  ${COMMON_LIBS}
  MyMPILib
  ${MPI_Fortran_LIBRARIES}
)

add_executable(test_binsrc.x source/test_binsrc.f90)
target_link_libraries(test_binsrc.x
  ${COMMON_LIBS}
)

add_executable(test_boozer_class.x source/test_boozer_class.f90)
target_link_libraries(test_boozer_class.x
  ${COMMON_LIBS}
)

add_executable(test_efit_class.x source/test_efit_class.f90)
target_link_libraries(test_efit_class.x
  ${COMMON_LIBS}
)

add_executable(test_geqdsk_tools.x source/test_geqdsk_tools.f90)
target_link_libraries(test_geqdsk_tools.x
  ${MAIN_LIB}
  ${MAGFIE_LIB}
)

if(LIBNEO_TESTING_ENABLED)
  set(TEST_OUTPUT_DIR ${CMAKE_BINARY_DIR}/test)
  file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR})
  # Create subdirectories for different test types
  file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR}/source)
  file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR}/scripts)
  file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR}/python)
  file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR}/poincare)

  add_executable(test_geoflux.x source/test_geoflux.f90)
  target_link_libraries(test_geoflux.x
    ${COMMON_LIBS}
    ${MAGFIE_LIB}
  )

  add_executable(geoflux_coord_dump.x source/geoflux_coord_dump.f90)
  target_link_libraries(geoflux_coord_dump.x
    ${COMMON_LIBS}
    ${MAGFIE_LIB}
  )
endif()

add_executable(test_simpson.x source/test_simpson.f90)
target_link_libraries(test_simpson.x
  ${COMMON_LIBS}
)
add_executable(test_hdf5_tools.x source/test_hdf5_tools.f90)
target_link_libraries(test_hdf5_tools.x
  ${HDF5_TOOLS_LIB}
)

add_executable(test_mympilib.x
  source/test_mympilib.f90
  source/derived_scheduler_module.f90
)
target_link_libraries(test_mympilib.x
  MyMPILib
  ${MPI_Fortran_LIBRARIES}
)

add_executable(test_system_utility.x
  source/test_system_utility.f90
  ../src/local_rusage.c)
target_link_libraries(test_system_utility.x
  ${COMMON_LIBS}
)

add_executable(test_util.x source/test_util.f90)
target_link_libraries(test_util.x ${COMMON_LIBS})

add_executable(test_interpolate.x source/test_interpolate.f90)
target_link_libraries(test_interpolate.x ${COMMON_LIBS})

add_executable(test_batch_interpolate.x source/test_batch_interpolate.f90)
target_link_libraries(test_batch_interpolate.x ${COMMON_LIBS})

add_executable(test_vmec_modules.x source/test_vmec_modules.f90)
target_link_libraries(test_vmec_modules.x neo)

add_executable(test_coordinate_systems.x source/test_coordinate_systems.f90)
target_link_libraries(test_coordinate_systems.x neo ${MAGFIE_LIB})

add_executable(test_chartmap_coordinates.x source/test_chartmap_coordinates.f90)
target_link_libraries(test_chartmap_coordinates.x neo ${MAGFIE_LIB})

add_executable(test_chartmap_validator.x source/test_chartmap_validator.f90)
target_link_libraries(test_chartmap_validator.x neo)

# TODO: Update test for new psi_1..psi_7 interface
# add_executable(test_analytical_gs_circular.x source/test_analytical_gs_circular.f90)
# target_link_libraries(test_analytical_gs_circular.x ${MAGFIE_LIB} ${COMMON_LIBS})

add_executable(test_analytical_circular.x source/test_analytical_circular.f90)
target_link_libraries(test_analytical_circular.x ${MAGFIE_LIB} ${COMMON_LIBS})

add_executable(test_nctools_nc_get3.x source/test_nctools_nc_get3.f90)
target_link_libraries(test_nctools_nc_get3.x
  ${COMMON_LIBS}
)


## Standard tests according to standard libneo Fortran test convention.

add_test(NAME test_binsrc COMMAND test_binsrc.x)
add_test(NAME test_boozer_class COMMAND test_boozer_class.x)
add_test(NAME test_efit_class COMMAND test_efit_class.x)
add_test(NAME test_geqdsk_tools COMMAND test_geqdsk_tools.x)
add_test(NAME test_simpson COMMAND test_simpson.x)
add_test(NAME test_hdf5_tools COMMAND test_hdf5_tools.x)
add_test(NAME test_util COMMAND test_util.x)
add_test(NAME test_interpolate COMMAND test_interpolate.x)
add_test(NAME test_batch_interpolate COMMAND test_batch_interpolate.x)
add_test(NAME test_collision_freqs COMMAND test_collision_freqs.x)
add_test(NAME test_transport COMMAND test_transport.x)
add_test(NAME test_vmec_modules COMMAND test_vmec_modules.x)
add_test(NAME test_coordinate_systems COMMAND test_coordinate_systems.x)
# add_test(NAME test_analytical_gs_circular COMMAND test_analytical_gs_circular.x)  # TODO: update for new interface
add_test(NAME test_analytical_circular COMMAND test_analytical_circular.x)
add_test(NAME test_nctools_nc_get3 COMMAND test_nctools_nc_get3.x)

add_test(
  NAME test_ascot5_compare
  COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/scripts/test_ascot5_compare.py
          --build-dir ${CMAKE_BINARY_DIR}
          --output-dir ${CMAKE_BINARY_DIR}/test/scripts
)
set_tests_properties(test_ascot5_compare PROPERTIES
  ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/test/scripts/ascot5/build:$ENV{LD_LIBRARY_PATH}"
)

if(LIBNEO_TESTING_ENABLED)
  add_test(
    NAME test_geoflux
    COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/scripts/run_geoflux_test.py
            --exe $<TARGET_FILE:test_geoflux.x>
            --data-dir ${TEST_OUTPUT_DIR}/scripts
  )

  add_test(
    NAME test_geqdsk_plot
    COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/scripts/geqdsk_plot.py
            --data-dir ${TEST_OUTPUT_DIR}/scripts
            --output-dir ${TEST_OUTPUT_DIR}/scripts
            --basename geqdsk
  )
  set_tests_properties(test_geqdsk_plot PROPERTIES DEPENDS test_geoflux)

  add_test(
    NAME test_geoflux_plot
    COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/scripts/geoflux_plot.py
            --exe $<TARGET_FILE:geoflux_coord_dump.x>
            --data-dir ${TEST_OUTPUT_DIR}/scripts
            --output-dir ${TEST_OUTPUT_DIR}/scripts
            --basename geoflux
  )
  set_tests_properties(test_geoflux_plot PROPERTIES DEPENDS test_geoflux)
endif()

## Custom tests

# This requires an additional call to setup_test_arnoldi.py, so use the script
# instead of calling the executable directly.
# \note Maybe this could be solved with setup_test_arnoldi.py && test_arnoldi.x.
# Afterwards fail condition is set: test will fail if the output
# contains 'STOP'.
add_test(NAME test_arnoldi_setup
  COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/source/setup_test_arnoldi.py)
set_tests_properties(test_arnoldi_setup PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_arnoldi COMMAND test_arnoldi.x)
set_tests_properties(test_arnoldi PROPERTIES 
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS test_arnoldi_setup
  FAIL_REGULAR_EXPRESSION "STOP")

add_test(NAME test_mympilib
  COMMAND test_mympilib.x)
set_tests_properties(test_mympilib PROPERTIES  PASS_REGULAR_EXPRESSION "Derived initMaster")
add_test(NAME test_system_utility
         COMMAND test_system_utility.x)
set_tests_properties(test_system_utility PROPERTIES  FAIL_REGULAR_EXPRESSION "WARNING: resource usage could not be determined.")

add_test(
  NAME setup_chartmap_volume
  COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/scripts/setup_chartmap_volume.py
          --output-dir ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(
  NAME setup_vmec_wout
  COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/scripts/setup_vmec_wout.py
          --output-dir ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(
  NAME setup_vmec_chartmap_python
  COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/scripts/setup_vmec_chartmap_python.py
          --output-dir ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(setup_vmec_chartmap_python PROPERTIES
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS setup_vmec_wout
)

add_test(
  NAME validate_vmec_map2disc_chartmap_boundary_python
  COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/scripts/validate_vmec_map2disc_chartmap_boundary.py
          --wout ${CMAKE_CURRENT_BINARY_DIR}/wout.nc
          --chartmap ${CMAKE_CURRENT_BINARY_DIR}/wout_vmec_python.chartmap.nc
)
set_tests_properties(validate_vmec_map2disc_chartmap_boundary_python PROPERTIES
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS setup_vmec_chartmap_python
)

add_test(
  NAME validate_vmec_map2disc_chartmap_boundary_cli
  COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/scripts/validate_vmec_map2disc_chartmap_boundary.py
          --wout ${CMAKE_CURRENT_BINARY_DIR}/wout.nc
          --chartmap ${CMAKE_CURRENT_BINARY_DIR}/wout_vmec_cli.chartmap.nc
)
set_tests_properties(validate_vmec_map2disc_chartmap_boundary_cli PROPERTIES
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS setup_vmec_chartmap_python
)

add_test(NAME test_chartmap_coordinates
  COMMAND test_chartmap_coordinates.x)
set_tests_properties(test_chartmap_coordinates PROPERTIES
  FAIL_REGULAR_EXPRESSION "STOP"
  DEPENDS setup_chartmap_volume)

add_test(NAME test_chartmap_validator
  COMMAND test_chartmap_validator.x)
set_tests_properties(test_chartmap_validator PROPERTIES
  FAIL_REGULAR_EXPRESSION "STOP"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_library(chartmap_test_utils STATIC source/chartmap_test_utils.f90)
target_link_libraries(chartmap_test_utils PUBLIC neo)

add_executable(test_vmec_chartmap_generator.x source/test_vmec_chartmap_generator.f90)
target_link_libraries(test_vmec_chartmap_generator.x chartmap_test_utils)
add_test(NAME test_vmec_chartmap_generator
  COMMAND test_vmec_chartmap_generator.x)
set_tests_properties(test_vmec_chartmap_generator PROPERTIES
  FAIL_REGULAR_EXPRESSION "STOP"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS setup_vmec_wout)

add_executable(test_python_chartmap_generator.x source/test_python_chartmap_generator.f90)
target_link_libraries(test_python_chartmap_generator.x chartmap_test_utils)
add_test(NAME test_python_chartmap_generator
  COMMAND test_python_chartmap_generator.x)
set_tests_properties(test_python_chartmap_generator PROPERTIES
  FAIL_REGULAR_EXPRESSION "STOP"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS setup_vmec_chartmap_python)

add_executable(test_vmec_coordinate_system.x source/test_vmec_coordinate_system.f90)
target_link_libraries(test_vmec_coordinate_system.x neo)
add_test(NAME test_vmec_coordinate_system
  COMMAND test_vmec_coordinate_system.x)
set_tests_properties(test_vmec_coordinate_system PROPERTIES
  FAIL_REGULAR_EXPRESSION "STOP"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS setup_vmec_wout)

add_executable(test_chartmap_matches_vmec.x source/test_chartmap_matches_vmec.f90)
target_link_libraries(test_chartmap_matches_vmec.x neo)
add_test(NAME test_chartmap_matches_vmec
  COMMAND test_chartmap_matches_vmec.x)
set_tests_properties(test_chartmap_matches_vmec PROPERTIES
  FAIL_REGULAR_EXPRESSION "STOP"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS setup_vmec_wout)

add_executable(test_chartmap_vmec_mapping.x source/test_chartmap_vmec_mapping.f90)
target_link_libraries(test_chartmap_vmec_mapping.x neo)
add_test(NAME test_chartmap_vmec_mapping
  COMMAND test_chartmap_vmec_mapping.x)
set_tests_properties(test_chartmap_vmec_mapping PROPERTIES
  FAIL_REGULAR_EXPRESSION "STOP"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS setup_vmec_wout)

add_executable(test_chartmap_derivatives.x source/test_chartmap_derivatives.f90)
target_link_libraries(test_chartmap_derivatives.x neo)
add_test(NAME test_chartmap_derivatives
  COMMAND test_chartmap_derivatives.x)
set_tests_properties(test_chartmap_derivatives PROPERTIES
  FAIL_REGULAR_EXPRESSION "STOP"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS setup_vmec_wout)

add_executable(test_refcoords_file_detection.x source/test_refcoords_file_detection.f90)
target_link_libraries(test_refcoords_file_detection.x neo)
add_test(NAME test_refcoords_file_detection
  COMMAND test_refcoords_file_detection.x)
set_tests_properties(test_refcoords_file_detection PROPERTIES
  FAIL_REGULAR_EXPRESSION "STOP"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS "setup_vmec_wout;setup_chartmap_volume")

add_subdirectory(util_for_test)
add_subdirectory(util_for_test_field)
add_subdirectory(field)
add_subdirectory(magfie)
add_subdirectory(polylag)
add_subdirectory(poincare)
add_subdirectory(odeint)
