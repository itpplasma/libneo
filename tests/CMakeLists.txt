
find_package(MPI REQUIRED)

set (MAIN_LIB ${MAIN_LIB} libneo)
set (HDF5_TOOLS_LIB ${HDF5_TOOLS_LIB} hdf5_tools)

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
find_package(MPI REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/src/MyMPILib)
include_directories(${CMAKE_SOURCE_DIR}/src/hdf5_tools)
include_directories(${MPI_Fortran_INCLUDE_PATH})


find_package(HDF5 COMPONENTS C Fortran Fortran_HL REQUIRED)
if (NOT ${HDF5_FOUND})
	message(FATAL_ERROR "HDF5 not found, please set environment variable HDF5_ROOT to the correct path")
endif ()
include_directories(${HDF5_INCLUDE_DIRS})

execute_process(COMMAND nc-config --flibs
                OUTPUT_VARIABLE NETCDF_LIBRARIES
                OUTPUT_STRIP_TRAILING_WHITESPACE)

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
  ${NETCDF_LIBRARIES}
  )

add_executable(test_hdf5_tools.x source/test_hdf5_tools.f90)
target_link_libraries(test_hdf5_tools.x
  ${HDF5_TOOLS_LIB}
  ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES}
  )

add_executable(test_mympilib.x
  source/test_mympilib.f90
  source/derived_scheduler_module.f90
  )
target_link_libraries(test_mympilib.x
  MyMPILib
  ${MPI_Fortran_LIBRARIES}
  )
