
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

find_package(BLAS)
find_package(LAPACK)
include_directories(${CMAKE_SOURCE_DIR}/src/MyMPILib)

add_executable(test_arnoldi.x source/test_arnoldi.f90)

target_link_libraries(test_arnoldi.x
  ${MAIN_LIB}
  MyMPILib
  ${MPI_Fortran_LIBRARIES}
  ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
  )

add_executable(test_efit_class.x source/test_efit_class.f90)
target_link_libraries(test_efit_class.x
  ${MAIN_LIB}
  ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
  )

add_executable(test_boozer_class.x source/test_boozer_class.f90)
target_link_libraries(test_boozer_class.x
  ${MAIN_LIB}
  ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
  )

add_executable(test_system_utility.x
  source/test_system_utility.f90
  ../src/local_rusage.c)
target_link_libraries(test_system_utility.x
  ${MAIN_LIB}
  ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
  )
