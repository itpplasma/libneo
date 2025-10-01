### Basic settings
cmake_minimum_required (VERSION 3.22)

### Define the project
# project(hdf5_tools)
enable_language(C Fortran)

### Set compiler flags
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-O2 -cpp -fPIC -Wall -fcheck=bounds -std=f2008")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-O2 -cpp -mkl -stand f08")
endif ()

### Find HDF5
find_package(HDF5 COMPONENTS C HL Fortran Fortran_HL REQUIRED)

add_library(hdf5_tools STATIC
  hdf5_tools.F90
  hdf5_tools_f2003.f90
)
target_link_libraries(hdf5_tools PRIVATE
  hdf5::hdf5
  hdf5::hdf5_hl
  hdf5::hdf5_fortran
  hdf5::hdf5_hl_fortran
)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  install (TARGETS hdf5_tools DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/lib)
  install (DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  install (TARGETS hdf5_tools DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ifort/lib)
  install (DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ifort/include)
endif ()

add_library(LIBNEO::hdf5_tools ALIAS hdf5_tools)
