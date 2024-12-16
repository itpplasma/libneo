### Basic settings
cmake_minimum_required(VERSION 3.18)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

### Define the project
project (hdf5_file_merging_tool)
enable_language(Fortran)

### Define the name of the executable
set(PROJECT_EXE_NAME "h5merge.x")

### Set compiler flags
set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)
set(CMAKE_Fortran_FLAGS "-g -Wall")

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(HDF5 COMPONENTS Fortran Fortran_HL)
if (NOT ${HDF5_FOUND})
	message(FATAL_ERROR "HDF5 not found, please set environment variable HDF5_ROOT to the correct path")
endif ()
include_directories(${HDF5_INCLUDE_DIRS})

find_library(HDF5_TOOLS_LIB NAMES hdf5_tools libhdf5_tools.a PATHS ${CMAKE_SOURCE_DIR}/../../build/ ${CMAKE_BINARY_DIR}/ PATH_SUFFIXES src/hdf5_tools/ NO_DEFAULT_PATH)

include_directories(${CMAKE_SOURCE_DIR}/../../build/include ${CMAKE_BINARY_DIR})

set(Interface_Libs
  ${HDF5_TOOLS_LIB}
  ${HDF5_Fortran_HL_LIBRARIES}
  ${HDF5_Fortran_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
)

### Define executable
add_executable(${PROJECT_EXE_NAME}
  ../src/h5merge.f90
)

### Linking
target_link_libraries(${PROJECT_EXE_NAME}
  ${Interface_Libs}
)

### For make install
# install (TARGETS ${PROJECT_EXE_NAME} DESTINATION /temp/)
# set_property(TARGET ${PROJECT_EXE_NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
