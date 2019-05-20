
find_package(MPI REQUIRED)

set (MAIN_LIB ${MAIN_LIB} libneo)

set(PROJLIBS /proj/plasma/Libs/ CACHE STRING "Common library path")
set(MYMPILIB_PATH ${PROJLIBS}/MyMPILib/Build-Release/)
include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories(${PROJECT_BINARY_DIR})
include_directories(${MYMPILIB_PATH}/OBJS)

### Path to OpenBLAS
set(OPEN_BLAS_PATH /usr/lib/openblas-base/)
set(OPEN_BLAS_LAPACK_LIB ${OPEN_BLAS_PATH})
set(OPEN_BLAS_LIB ${OPEN_BLAS_PATH})


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