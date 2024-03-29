### Basic settings
cmake_minimum_required (VERSION 3.12)

### Define the project
# project(hdf5_tools)
enable_language(Fortran)

### Set compiler flags
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-O2 -cpp -fPIC -Wall -fcheck=bounds -std=f2008")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-O2 -cpp -mkl -stand f08")
endif ()
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

### Find HDF5
find_package(HDF5 COMPONENTS Fortran Fortran_HL)
if (NOT ${HDF5_FOUND})
	message(FATAL_ERROR "HDF5 not found, please set environment variable HDF5_ROOT to the correct path")
endif ()

# Set variables in parent scope, so find_package is not needed there
set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS} PARENT_SCOPE)
set(HDF5_LIBRARIES ${HDF5_LIBRARIES} PARENT_SCOPE)
set(HDF5_HL_LIBRARIES ${HDF5_HL_LIBRARIES} PARENT_SCOPE)


include_directories(${HDF5_INCLUDE_DIRS})

### Define executable
add_library(hdf5_tools
	${CMAKE_SOURCE_DIR}/src/hdf5_tools/hdf5_tools.f90
	${CMAKE_SOURCE_DIR}/src/hdf5_tools/hdf5_tools_f2003.f90
)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  install (TARGETS hdf5_tools DESTINATION ${CMAKE_SOURCE_DIR}/lib)
  install (DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${CMAKE_SOURCE_DIR}/include)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  install (TARGETS hdf5_tools DESTINATION ${CMAKE_SOURCE_DIR}/ifort/lib)
  install (DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${CMAKE_SOURCE_DIR}/ifort/include)
endif ()
