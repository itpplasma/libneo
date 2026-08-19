cmake_minimum_required(VERSION 3.22)
project(hdf5_file_merging_tool LANGUAGES C Fortran)

include(FetchContent)
set(BUILD_TESTING OFF CACHE BOOL "Do not build dependency tests" FORCE)
FetchContent_Declare(
  fortio
  GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
  GIT_TAG 471256b7382afb8f641f37d833cc9172605da018
)
FetchContent_MakeAvailable(fortio)

add_executable(h5merge.x src/h5merge.f90)
target_link_libraries(h5merge.x PRIVATE fortio::fortio)
