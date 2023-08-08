cmake_minimum_required(VERSION 3.12)
# project (MyMPILib)
enable_language (Fortran)

# MPI
find_package(MPI REQUIRED)

include (${CMAKE_SOURCE_DIR}/src/MyMPILib/ProjectConfig.cmake.in)
option(MPE_SUPPORT "Should the library be built with MPE-Support?" OFF)

message(STATUS "The Compiler ID is ${CMAKE_Fortran_COMPILER_ID}")

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

	#gfortran -I/usr/lib/openmpi/include -pthread -I/usr/lib/openmpi/lib -L/usr//lib -L/usr/lib/openmpi/lib -lmpi_f90 -lmpi_f77 -lmpi -ldl -lhwloc

#	include_directories(${INCLUDEDIRS_MPI_GNU})
#	link_directories(${LINKDIRS_MPI_GNU})

	set (CMAKE_Fortran_FLAGS "-cpp -pthread")
        #set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
        #set (CMAKE_Fortran_FLAGS_DEBUG "")

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

	#ifort -I/usr/local/openmpi-1.6.3-intel/include -I/usr/local/openmpi-1.6.3-intel/lib -L/usr/local/openmpi-1.6.3-intel/lib -lmpi_f90 -lmpi_f77 -lmpi -ldl -lm -lnuma -Wl,--export-dynamic -lrt -lnsl -lutil

#	include_directories(${INCLUDEDIRS_MPI_INTEL})
#	link_directories(${LINKDIRS_MPI_INTEL})

	set (CMAKE_Fortran_FLAGS "-cpp -heap-arrays")
        #set (CMAKE_Fortran_FLAGS_RELEASE "")
        #set (CMAKE_Fortran_FLAGS_DEBUG "")
endif ()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${CMAKE_SOURCE_DIR}/src/MyMPILib/Generic")
include_directories("${CMAKE_SOURCE_DIR}/src/MyMPILib/Specific")
include_directories("${CMAKE_SOURCE_DIR}/src/MyMPILib/Internal")
include_directories("${CMAKE_SOURCE_DIR}/src/MyMPILib/Tools")

if (MPE_SUPPORT)
	message(STATUS "MPE-Support activated, MPE-Path: ${MPE_PATH}")
	if (DEFINED MPE_PATH)
		include_directories(${MPE_PATH}/include)
	endif ()
	add_definitions(-DMPE_SUPPORT)
endif ()

add_library(MyMPILib
${CMAKE_SOURCE_DIR}/src/MyMPILib/Specific/packBuffer_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Specific/mpiprovider_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Tools/list_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Internal/packable_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Internal/matrix_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Internal/clientStatus_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Generic/workunit_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Internal/wuDataRequester_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Tools/intList_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Tools/wuList_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Tools/packableList_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Tools/commandline_parser_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Generic/scheduler_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Tools/mpelog_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Tools/myLog_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Tools/configFile_parser_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Internal/wuMergeWorkunit_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Internal/wuMergeChunk_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Generic/initWorkunit_module.f90
${CMAKE_SOURCE_DIR}/src/MyMPILib/Generic/genericWorkunit_module.f90
)

install(TARGETS MyMPILib ARCHIVE DESTINATION lib)

add_custom_target(
    Versioning
    COMMAND ${CMAKE_SOURCE_DIR}/src/MyMPILib/Scripts/do_versioning.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/MyMPILib
)

add_dependencies(MyMPILib Versioning)
