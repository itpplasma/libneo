# This file is intended for configuration of libneo.
# List of source files are in CMakeSources.in.

# Basic cmake settings/requirements
cmake_minimum_required(VERSION 3.20)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Basic project settings.
project(libneo)

# The version number.
set (libneo_VERSION_MAJOR 0)
set (libneo_VERSION_MINOR 1)

enable_language(Fortran)
enable_testing()

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# set (CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

if (EXISTS ${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
  include (${CMAKE_BINARY_DIR}/BuildConfig.cmake.in)
endif()

find_program (BASH_PROGRAM bash)

find_package (MPI REQUIRED)

# Include file for the sources.
include(CMakeSources.in)


if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

  set (CMAKE_Fortran_FLAGS "-cpp -g -Og -fcheck=all -fbacktrace \
    -fno-realloc-lhs -fopenmp -std=f2008 -Wall -Wextra")

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

  set (CMAKE_Fortran_FLAGS "-cpp -check all -traceback \
    -assume norealloc_lhs -stand f08 -qopenmp -warn all")

endif ()

include_directories(${MPI_Fortran_INCLUDE_PATH})
message(STATUS "The MPI_Fortran_INCLUDE_PATH is ${MPI_Fortran_INCLUDE_PATH}")

find_program(NF_CONFIG "nf-config")

if (NF_CONFIG)
execute_process(COMMAND nf-config --includedir
                OUTPUT_VARIABLE NETCDFINCLUDE_DIR)
execute_process(COMMAND nc-config --libdir
				OUTPUT_VARIABLE NETCDFLIB_DIR)
execute_process(COMMAND nf-config --flibs
                OUTPUT_VARIABLE NETCDF_FLIBS)
else()
message(SEND_ERROR "nf-config not found. Please install libnetcdff-dev")
endif()

string(STRIP ${NETCDFINCLUDE_DIR} NETCDFINCLUDE_DIR)
string(STRIP ${NETCDFLIB_DIR} NETCDFLIB_DIR)
string(STRIP ${NETCDF_FLIBS} NETCDF_FLIBS)

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

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

#set_target_properties(libneo PROPERTIES VERSION ${libneo_VERSION})

add_subdirectory(src)
add_subdirectory(tests)

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

# Python bindings
cmake_minimum_required(VERSION 3.17.2...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)

# F2PY headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# Print F2PY headers
message(STATUS "F2PY headers: ${F2PY_INCLUDE_DIR}")

#add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
#target_link_libraries(fortranobject PUBLIC Python::NumPy)
#target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
#set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

#add_custom_command(
#  OUTPUT interpolatemodule.c interpolate-f2pywrappers.f90
#  DEPENDS src/interpolate.f90
#  VERBATIM
#  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
#  "${CMAKE_CURRENT_SOURCE_DIR}/src/interpolate.f90" -m example --lower)

#python_add_library(example MODULE "${CMAKE_CURRENT_BINARY_DIR}/interpolatemodule.c"
#                   "${CMAKE_CURRENT_SOURCE_DIR}/src/interpolate.f90" WITH_SOABI)
#target_link_libraries(example PRIVATE fortranobject)

#install(TARGETS example DESTINATION .)
