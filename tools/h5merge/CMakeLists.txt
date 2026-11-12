cmake_minimum_required(VERSION 3.22)
project(hdf5_file_merging_tool LANGUAGES C Fortran)

include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/ItpDeps.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/DependencyPins.cmake)
itp_add_dependency(fortio REF ${LIBNEO_FORTIO_REF}
  GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
  TARGET fortio)

add_executable(h5merge.x src/h5merge.f90)
target_link_libraries(h5merge.x PRIVATE fortio::fortio)
