set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

find_package(MPI REQUIRED)

set (MAIN_LIB ${MAIN_LIB} libneo)
set (HDF5_TOOLS_LIB ${HDF5_TOOLS_LIB} hdf5_tools)
set (MAGFIE_LIB ${MAGFIE_LIB} magfie)

set(PROJLIBS /proj/plasma/Libs/ CACHE STRING "Common library path")
include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories(${PROJECT_BINARY_DIR})


if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

  set (CMAKE_Fortran_FLAGS "-g -cpp -std=f2008")

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

  set (CMAKE_Fortran_FLAGS "-g -cpp -heap-arrays -stand f08")

endif ()

find_package(BLAS)
find_package(LAPACK)
find_package(MPI REQUIRED)

include_directories(${CMAKE_BINARY_DIR}/src/)
include_directories(${CMAKE_SOURCE_DIR}/src/MyMPILib)
include_directories(${CMAKE_SOURCE_DIR}/src/hdf5_tools)
include_directories(${MPI_Fortran_INCLUDE_PATH})


find_package(HDF5 COMPONENTS C Fortran Fortran_HL REQUIRED)
if (NOT ${HDF5_FOUND})
	message(FATAL_ERROR "HDF5 not found, please set environment variable HDF5_ROOT to the correct path")
endif ()
include_directories(${HDF5_INCLUDE_DIRS})

execute_process(COMMAND nc-config --flibs
                OUTPUT_VARIABLE NETCDF_LIBRARIES
                OUTPUT_STRIP_TRAILING_WHITESPACE)

set (COMMON_LIBS ${MAIN_LIB} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES})

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

add_executable(test_hdf5_tools.x source/test_hdf5_tools.f90)
target_link_libraries(test_hdf5_tools.x
  ${HDF5_TOOLS_LIB}
  ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES}
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


## Standard tests according to standard libneo Fortran test convention.

function (add_fortran_test name)
  add_test(NAME ${name} COMMAND ${name}.x)
  set_tests_properties(${name} PROPERTIES  FAIL_REGULAR_EXPRESSION "STOP")
endfunction()

add_fortran_test(test_binsrc)
add_fortran_test(test_boozer_class)
add_fortran_test(test_efit_class)
add_fortran_test(test_geqdsk_tools)
add_fortran_test(test_hdf5_tools)
add_fortran_test(test_util)

## Custom tests

# This requires an additional call to test_arnoldi.py, so use the script
# instead of calling the executable directly.
# \note Maybe this could be solved with test_arnoldi.py && test_arnoldi.x.
# Afterwards fail condition is set: test will fail if the output
# contains 'STOP'.
add_test(NAME test_arnoldi
         COMMAND ../../tests/python_scripts/test_arnoldi.py && test_arnoldi.x)
set_tests_properties(test_arnoldi PROPERTIES  FAIL_REGULAR_EXPRESSION "STOP")

add_test(NAME test_mympilib
  COMMAND test_mympilib.x)
set_tests_properties(test_mympilib PROPERTIES  PASS_REGULAR_EXPRESSION "Derived initMaster")
add_test(NAME test_system_utility
         COMMAND test_system_utility.x)
set_tests_properties(test_system_utility PROPERTIES  FAIL_REGULAR_EXPRESSION "WARNING: resource usage could not be determined.")
