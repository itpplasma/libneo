
find_package(MPI REQUIRED)

set (MAIN_LIB ${MAIN_LIB} libneo)

set(PROJLIBS /proj/plasma/Libs/ CACHE STRING "Common library path")
include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories(${PROJECT_BINARY_DIR})


if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

  set (CMAKE_Fortran_FLAGS "-cpp -std=f2008")

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

  set (CMAKE_Fortran_FLAGS "-cpp -heap-arrays -stand f08")

endif ()

# Settings based on OS version.
### Path to OpenBLAS
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL 8.3.0)
    set(MYMPILIB_PATH ${PROJLIBS}/MyMPILib/Build-gnu-8.3.0/)
    set(OPEN_BLAS_PATH /usr/lib/openblas/ /usr/lib/x86_64-linux-gnu/)
  else()
    set(MYMPILIB_PATH ${PROJLIBS}/MyMPILib/Build-Release/)
    set(OPEN_BLAS_PATH /usr/lib/openblas-base/ /usr/lib/x86_64-linux-gnu/)
  endif()
else()
  set(MYMPILIB_PATH ${PROJLIBS}/MyMPILib/Build-Release/)
  set(OPEN_BLAS_PATH /usr/lib/openblas-base/ /usr/lib/x86_64-linux-gnu/)
endif()
set(OPEN_BLAS_LAPACK_LIB ${OPEN_BLAS_PATH})
set(OPEN_BLAS_LIB ${OPEN_BLAS_PATH})

include_directories(${MYMPILIB_PATH}/OBJS)

find_library(open_blas_lib blas ${OPEN_BLAS_LIB} NO_DEFAULT_PATH)
find_library(open_blas_lapack_lib lapack ${OPEN_BLAS_LAPACK_LIB} NO_DEFAULT_PATH)
find_library(MyMPILib_lib MyMPILib ${MYMPILIB_PATH})

set(LINALG_LIBRARIES ${open_blas_lib} ${open_blas_lapack_lib})

add_executable(test_arnoldi.x source/test_arnoldi.f90)

target_link_libraries(test_arnoldi.x
  ${MAIN_LIB}
  ${MyMPILib_lib}
  ${MPI_Fortran_LIBRARIES}
  ${LINALG_LIBRARIES}
  )

add_executable(test_efit_class.x source/test_efit_class.f90)
target_link_libraries(test_efit_class.x
  ${MAIN_LIB}
  ${LINALG_LIBRARIES}
  )

add_executable(test_boozer_class.x source/test_boozer_class.f90)
target_link_libraries(test_boozer_class.x
  ${MAIN_LIB}
  ${LINALG_LIBRARIES}
  )

add_executable(test_system_utility.x
  source/test_system_utility.f90
  ../src/local_rusage.c)
target_link_libraries(test_system_utility.x
  ${MAIN_LIB}
  ${LINALG_LIBRARIES}
  )
