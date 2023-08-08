# This file is intended for configuration of libneo.
# List of source files are in CMakeSources.in.

# Basic cmake settings/requirements
cmake_minimum_required(VERSION 3.12)
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

find_program(NC_CONFIG "nf-config")

if (NC_CONFIG)
execute_process(COMMAND nf-config --includedir
                OUTPUT_VARIABLE NETCDFINCLUDE_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND nc-config --libdir
                OUTPUT_VARIABLE NETCDFLIB_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND nf-config --flibs
                OUTPUT_VARIABLE NETCDF_FLIBS
                OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
message(SEND_ERROR "nc-config not found. Please install libnetcdff-dev")
endif()

message(STATUS "NetCDF include path: " ${NETCDFINCLUDE_DIR})
message(STATUS "NetCDF lib path: " ${NETCDFLIB_DIR})
message(STATUS "NetCDF Fortran libs: " ${NETCDF_FLIBS})

# Replace space by semicolon in the Fortran libs
string(REPLACE " " ";" NETCDF_FLIBS ${NETCDF_FLIBS})

include_directories(${NETCDFINCLUDE_DIR})
link_directories(${NETCDFLIB_DIR})
add_link_options(${NETCDF_FLIBS})

add_subdirectory(src/MyMPILib)
add_subdirectory(src/hdf5_tools)

# RKF45 from https://people.sc.fsu.edu/~jburkardt/
# licenced under LGPL (see SRC/contrib/LICENSE.rkf45)
add_library ( rkf45 SHARED
  src/contrib/rkf45.f90
)

# Main goal (of the created build file): creating the library file.
add_library(libneo SHARED ${SOURCE_FILES})

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(libneo
  rkf45 ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

#set_target_properties(libneo PROPERTIES VERSION ${libneo_VERSION})

add_subdirectory(src)
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
         COMMAND "tests/test_hdf5_tools.x")
set_tests_properties(Test_hdf5_tools PROPERTIES  FAIL_REGULAR_EXPRESSION "STOP")

add_test(NAME Test_mympilib
         COMMAND tests/test_mympilib.x)
set_tests_properties(Test_mympilib PROPERTIES  PASS_REGULAR_EXPRESSION "Derived initMaster")
add_test(NAME Test_system_utility
         COMMAND tests/test_system_utility.x)
set_tests_properties(Test_system_utility PROPERTIES  FAIL_REGULAR_EXPRESSION "WARNING: resource usage could not be determined.")

# MacOS RPATH specifics
# see https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")


#~ # For testing purposes, e.g. printing the environment.
#~ add_test(NAME run_xterm
         #~ COMMAND xterm)
