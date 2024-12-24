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

### Find HDF5
find_package(HDF5 COMPONENTS Fortran Fortran_HL)
if (NOT ${HDF5_FOUND})
	message(FATAL_ERROR "HDF5 not found, please set environment variable HDF5_ROOT to the correct path")
endif ()

if (APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
  message(STATUS "Apple M1 detected, fixing paths to Fortran HDF5 libraries")
  set(HDF5_LIBRARIES /opt/homebrew/lib/libhdf5_fortran.dylib)
  set(HDF5_HL_LIBRARIES /opt/homebrew/lib/libhdf5_hl_fortran.dylib)
endif()

# Set variables in parent scope, so find_package is not needed there
set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS} PARENT_SCOPE)
set(HDF5_LIBRARIES ${HDF5_LIBRARIES} PARENT_SCOPE)
set(HDF5_HL_LIBRARIES ${HDF5_HL_LIBRARIES} PARENT_SCOPE)


message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")
message(STATUS "HDF5 HL libraries: ${HDF5_HL_LIBRARIES}")


include_directories(${HDF5_INCLUDE_DIRS})

### Define executable
add_library(hdf5_tools
	        hdf5_tools.f90
	        hdf5_tools_f2003.f90
)
target_link_libraries(hdf5_tools PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  install (TARGETS hdf5_tools DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/lib)
  install (DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  install (TARGETS hdf5_tools DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ifort/lib)
  install (DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ifort/include)
endif ()

add_library(LIBNEO::hdf5_tools ALIAS hdf5_tools)
