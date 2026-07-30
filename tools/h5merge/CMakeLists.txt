cmake_minimum_required(VERSION 3.22)
project(hdf5_file_merging_tool LANGUAGES C Fortran)

include(FetchContent)
set(BUILD_TESTING OFF CACHE BOOL "Do not build dependency tests" FORCE)
FetchContent_Declare(
  fortio
  GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
  GIT_TAG bca5cd8234d17446aea09526bfb45736f96b5746
)
FetchContent_MakeAvailable(fortio)

add_executable(h5merge.x src/h5merge.f90)
target_link_libraries(h5merge.x PRIVATE fortio::fortio)
