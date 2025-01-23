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
find_package(HDF5 COMPONENTS C HL Fortran Fortran_HL)

include_directories(${CMAKE_SOURCE_DIR}/../../build/include ${CMAKE_BINARY_DIR})

set(Interface_Libs
  hdf5::hdf5
  hdf5::hdf5_hl
  hdf5::hdf5_fortran
  hdf5::hdf5_hl_fortran
  LAPACK::LAPACK
  BLAS::BLAS
)

### Define executable
add_executable(${PROJECT_EXE_NAME}
  ../src/h5merge.f90
)

### Linking
target_link_libraries(${PROJECT_EXE_NAME}
  ${CMAKE_SOURCE_DIR}/../../build/libneo.a
  ${CMAKE_SOURCE_DIR}/../../build/src/hdf5_tools/libhdf5_tools.a
  ${Interface_Libs}
)


### For make install
# install (TARGETS ${PROJECT_EXE_NAME} DESTINATION /temp/)
# set_property(TARGET ${PROJECT_EXE_NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
