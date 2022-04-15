# This file is intended for configuration of libneo.
# List of source files are in CMakeSources.in.

# Basic cmake settings/requirements
cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Basic project settings.
project(libneo)

# The version number.
set (libneo_VERSION_MAJOR 0)
set (libneo_VERSION_MINOR 1)

enable_language(Fortran)
enable_testing()

# set (CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

if (EXISTS ${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
  include (${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
endif()

find_program (BASH_PROGRAM bash)

find_package (MPI REQUIRED)

# Include file for the sources.
include(CMakeSources.in)


if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

  set (CMAKE_Fortran_FLAGS "-cpp -g -fbacktrace -std=f2008")

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

  set (CMAKE_Fortran_FLAGS "-cpp -heap-arrays -stand f08")

endif ()

include_directories(${MPI_Fortran_INCLUDE_PATH})
message(STATUS "The MPI_Fortran_INCLUDE_PATH is ${MPI_Fortran_INCLUDE_PATH}")
execute_process(COMMAND nc-config --includedir
                OUTPUT_VARIABLE NETCDFINCLUDE_DIR)
include_directories(${NETCDFINCLUDE_DIR})

add_subdirectory(src/MyMPILib)
add_subdirectory(src/hdf5_tools)

# Main goal (of the created build file): creating the library file.
add_library(libneo SHARED ${SOURCE_FILES})
#target_link_libraries(libneo ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_library(mc_efit SHARED
  src/magfie/polylag_3.f
  src/magfie/polylag_5.f
  src/magfie/spline5_RZ.f90
  src/magfie/amn_mod.f90
  src/magfie/bdivfree.f90
  src/magfie/bdivfree_mod.f90
  src/magfie/extract_fluxcoord_mod.f90
  src/magfie/field_c_mod.f90
  src/magfie/field_divB0.f90
  src/magfie/field_eq_mod.f90
  src/magfie/field_mesh3d.f90
  src/magfie/field_mod.f90
  src/magfie/input_files.f90
  src/magfie/inthecore_mod.f90
  src/magfie/magfie_cyl.f90
  src/magfie/magfield_mod.f90
  src/magfie/neo_magfie_mod.f90
  src/magfie/period_mod.f90
  src/magfie/theta_rz_mod.f90
  src/contrib/rkf45.f90
  src/odeint_rkf45.f90
  src/MC/alpha_lifetime_mod.f90
  src/MC/chamb_divB0.f90
  src/MC/rng.f90
  src/MC/collis_nbi.f90
  src/MC/sub_alpha_lifetime.f90)

#set_target_properties(libneo PROPERTIES VERSION ${libneo_VERSION})

add_subdirectory(tests)

add_test(Test_Build_Code make)
# This requires an additional call to test_arnoldi.py, so use the script
# instead of calling the executable directly.
# \note Maybe this could be solved with test_arnoldi.py && test_arnoldi.x.
# Afterwards fail condition is set: test will fail if the output
# contains 'STOP'.
add_test(NAME Test_arnoldi
         COMMAND ../tests/python_scripts/test_arnoldi.py && tests/test_arnoldi.x)
set_tests_properties(Test_arnoldi PROPERTIES  FAIL_REGULAR_EXPRESSION "STOP")
add_test(NAME Test_boozer
         COMMAND tests/test_boozer_class.x)
set_tests_properties(Test_boozer PROPERTIES  FAIL_REGULAR_EXPRESSION "STOP")
add_test(NAME Test_efit
         COMMAND tests/test_efit_class.x)
set_tests_properties(Test_efit PROPERTIES  FAIL_REGULAR_EXPRESSION "STOP")
add_test(NAME Test_hdf5_tools
         COMMAND "tests/test_hdf5_tools.x" && "rm test_file.h5")
set_tests_properties(Test_hdf5_tools PROPERTIES  FAIL_REGULAR_EXPRESSION "STOP")
add_test(NAME Test_mympilib
         COMMAND tests/test_mympilib.x)
set_tests_properties(Test_mympilib PROPERTIES  PASS_REGULAR_EXPRESSION "Derived initMaster")
add_test(NAME Test_system_utility
         COMMAND tests/test_system_utility.x)
set_tests_properties(Test_system_utility PROPERTIES  FAIL_REGULAR_EXPRESSION "WARNING: resource usage could not be determined.")

#~ # For testing purposes, e.g. printing the environment.
#~ add_test(NAME run_xterm
         #~ COMMAND xterm)
