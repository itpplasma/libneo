# This file is intended for configuration of libneo.
# List of source files are in CMakeSources.in.

# Basic cmake settings/requirements
cmake_minimum_required(VERSION 2.6)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Basic project settings.
project(libneo)

# The version number.
set (libneo_VERSION_MAJOR 0)
set (libneo_VERSION_MINOR 1)

enable_language(Fortran)
enable_testing()

#~ set (CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

if (EXISTS ${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
  include (${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
endif()

find_program (BASH_PROGRAM bash)

find_package (MPI REQUIRED)

# Include file for the sources.
include(CMakeSources.in)

set(PROJLIBS /proj/plasma/Libs/ CACHE STRING "Common library path")
set(MYMPILIB_PATH ${PROJLIBS}/MyMPILib/Build-Release/)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

  set (CMAKE_Fortran_FLAGS "-cpp -std=f2008")

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

  set (CMAKE_Fortran_FLAGS "-cpp -heap-arrays -stand f08")

endif ()

include_directories(${MPI_Fortran_INCLUDE_PATH})
message(STATUS "The MPI_Fortran_INCLUDE_PATH is ${MPI_Fortran_INCLUDE_PATH}")
include_directories(${MYMPILIB_PATH}/OBJS)

find_library(MyMPILib_lib MyMPILib ${MYMPILIB_PATH})

# Main goal (of the created build file): creating the library file.
add_library(libneo SHARED ${SOURCE_FILES})

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
#~ # For testing purposes, e.g. printing the environment.
#~ add_test(NAME run_xterm
         #~ COMMAND xterm)
