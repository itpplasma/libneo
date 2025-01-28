cmake_minimum_required(VERSION 3.12)

project (MyMPILib)
enable_language (C Fortran)

# MPI
find_package (MPI REQUIRED COMPONENTS Fortran)

include(${PROJECT_SOURCE_DIR}/ProjectConfig.cmake.in)
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

	set (CMAKE_Fortran_FLAGS "-cpp")
        #set (CMAKE_Fortran_FLAGS_RELEASE "")
        #set (CMAKE_Fortran_FLAGS_DEBUG "")
endif ()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/Generic")
include_directories("${PROJECT_SOURCE_DIR}/Specific")
include_directories("${PROJECT_SOURCE_DIR}/Internal")
include_directories("${PROJECT_SOURCE_DIR}/Tools")

if (MPE_SUPPORT)
	message(STATUS "MPE-Support activated, MPE-Path: ${MPE_PATH}")
	if (DEFINED MPE_PATH)
		include_directories(${MPE_PATH}/include)
	endif ()
	add_definitions(-DMPE_SUPPORT)
endif ()

add_library(MyMPILib
${PROJECT_SOURCE_DIR}/Specific/packBuffer_module.f90
${PROJECT_SOURCE_DIR}/Specific/mpiprovider_module.f90
${PROJECT_SOURCE_DIR}/Tools/list_module.f90
${PROJECT_SOURCE_DIR}/Internal/packable_module.f90
${PROJECT_SOURCE_DIR}/Internal/matrix_module.f90
${PROJECT_SOURCE_DIR}/Internal/clientStatus_module.f90
${PROJECT_SOURCE_DIR}/Generic/workunit_module.f90
${PROJECT_SOURCE_DIR}/Internal/wuDataRequester_module.f90
${PROJECT_SOURCE_DIR}/Tools/intList_module.f90
${PROJECT_SOURCE_DIR}/Tools/wuList_module.f90
${PROJECT_SOURCE_DIR}/Tools/packableList_module.f90
${PROJECT_SOURCE_DIR}/Tools/commandline_parser_module.f90
${PROJECT_SOURCE_DIR}/Generic/scheduler_module.f90
${PROJECT_SOURCE_DIR}/Tools/mpelog_module.f90
${PROJECT_SOURCE_DIR}/Tools/myLog_module.f90
${PROJECT_SOURCE_DIR}/Tools/configFile_parser_module.f90
${PROJECT_SOURCE_DIR}/Internal/wuMergeWorkunit_module.f90
${PROJECT_SOURCE_DIR}/Internal/wuMergeChunk_module.f90
${PROJECT_SOURCE_DIR}/Generic/initWorkunit_module.f90
${PROJECT_SOURCE_DIR}/Generic/genericWorkunit_module.f90
)
target_link_libraries(MyMPILib PRIVATE MPI::MPI_Fortran)
install(TARGETS MyMPILib ARCHIVE DESTINATION lib)

add_custom_target(
    Versioning
    COMMAND ${PROJECT_SOURCE_DIR}/Scripts/do_versioning.sh
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)

add_dependencies(MyMPILib Versioning)

add_library(LIBNEO::MyMPILib ALIAS MyMPILib)
