# This file is intended for configuration of libneo.
# List of source files are in CMakeSources.in.

# Basic cmake settings/requirements
cmake_minimum_required(VERSION 3.22)
if(APPLE AND NOT DEFINED HDF5_DIR)
  # locate brew
  find_program(BREW_EXECUTABLE brew HINTS /opt/homebrew/bin /usr/local/bin)
  if(BREW_EXECUTABLE)
    execute_process(
      COMMAND "${BREW_EXECUTABLE}" --prefix hdf5
      OUTPUT_VARIABLE _hdf5_prefix
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _hdf5_res
    )
    if(_hdf5_res EQUAL 0 AND EXISTS "${_hdf5_prefix}/lib/cmake/hdf5")
      set(HDF5_DIR "${_hdf5_prefix}/lib/cmake/hdf5" CACHE PATH "Path to HDF5Config.cmake")
      message(STATUS "HDF5_DIR set to ${HDF5_DIR}")
    endif()
  endif()
endif()

# Basic project settings.
project(libneo)

# Detect if libneo is standalone or a subproject
set(LIBNEO_IS_STANDALONE OFF)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(LIBNEO_IS_STANDALONE ON)
endif()

option(LIBNEO_BUILD_TESTING "Build tests for libneo" ${LIBNEO_IS_STANDALONE})

# Allow BUILD_TESTING to override (for subproject compatibility)
if(DEFINED BUILD_TESTING)
    set(LIBNEO_BUILD_TESTING ${BUILD_TESTING})
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules/")

# Define available build types
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Debug Release Fast
)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build." FORCE)
endif()

message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")

# Enable interprocedural optimization (LTO) for Release builds
# DISABLED: LTO causes numerical instability in geqdsk_GS_prefactor calculation
# This is a confirmed gfortran LTO bug affecting finite difference calculations
# set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)

include(Util)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# The version number.
set(libneo_VERSION_MAJOR 2025)
set(libneo_VERSION_MINOR 08)
set(libneo_VERSION_PATCH 20)

# Auto-detect CUDA for NVHPC compilers (must be before enable_language)
# This sets NVHPC_CUDA_HOME and NVCOMPILER_COMM_LIBS_HOME environment variables
# if NVHPC compiler is detected and CUDA is found at /opt/cuda or /usr/local/cuda
include(${PROJECT_SOURCE_DIR}/cmake/SetupNVHPCCuda.cmake)

enable_language(C Fortran)
if(LIBNEO_BUILD_TESTING)
    include(CTest)
endif()

option(ENABLE_OPENMP "Enable OpenMP compiler flags." ON)
option(ENABLE_OPENACC "Enable OpenACC compiler flags." ON)
set(OPENACC_OFFLOAD_TARGET "none" CACHE STRING "OpenACC offload target for GNU (none|nvptx)")
set_property(CACHE OPENACC_OFFLOAD_TARGET PROPERTY STRINGS none nvptx)
option(LIBNEO_ENABLE_GOLDEN_TESTS
    "Enable fetching and running the ODEINT golden record regression test" ON)

# Dependency fetching options
option(PREFER_SYSTEM_LIBS "Try system libraries before building from source" ON)
option(FORCE_FETCH_DEPS "Always build Fortran dependencies from source" OFF)

if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

find_program(BASH_PROGRAM bash)

find_package(MPI REQUIRED COMPONENTS C Fortran)
# Compiler-specific base flags and build profiles
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    # Base flags for all builds
    add_compile_options(-fPIC -g)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-cpp>)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fno-realloc-lhs>)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fmax-errors=1>)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fbacktrace>)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-132>)
    if(ENABLE_OPENACC)
        add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fopenacc>)
        add_link_options($<$<LINK_LANGUAGE:Fortran>:-fopenacc>)
        if(OPENACC_OFFLOAD_TARGET STREQUAL "nvptx")
            add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-foffload=nvptx-none>)
            add_link_options($<$<LINK_LANGUAGE:Fortran>:-foffload=nvptx-none>)
        elseif(NOT OPENACC_OFFLOAD_TARGET STREQUAL "none")
            message(FATAL_ERROR "Invalid OPENACC_OFFLOAD_TARGET='${OPENACC_OFFLOAD_TARGET}' (expected none|nvptx)")
        endif()
    endif()

    # Build type specific flags
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_compile_options(-Og
                          -Wall -Wextra
                          -Wno-unused-function)
        add_compile_options($<$<COMPILE_LANGUAGE:C>:-Wno-error=cast-function-type>)
        add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-fcheck=all>)
        add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Wno-external-argument-mismatch>)
        add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ffpe-trap=invalid,zero,overflow>)
    elseif(CMAKE_BUILD_TYPE MATCHES "Release|Fast")
        # Common optimization flags for Release and Fast builds
        add_compile_options(-O3 -DNDEBUG
                          -ffast-math -ffp-contract=fast
                          -funroll-loops
                          -ftree-vectorize)

        # Architecture-specific flags
        if(CMAKE_BUILD_TYPE MATCHES Fast OR CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
            # Fast mode or ARM: use native architecture (no portable baseline for ARM)
            add_compile_options(-march=native)
        else()
            # Release on x86-64: use v2 baseline for broad compatibility
            add_compile_options(-march=x86-64-v2 -mtune=generic)
        endif()
    endif()

    # Coverage disabled/removed from build system

    # https://github.com/numpy/numpy/issues/25777
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-error=incompatible-pointer-types")

    # Make trampolines a hard error - they would create executable stack
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Wtrampolines>)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Werror=trampolines>)
    # Apple's linker doesn't support -z,noexecstack (uses different security model)
    if(NOT APPLE)
        add_link_options(-Wl,-z,noexecstack)
    endif()

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    # Base flags for all builds
    add_compile_options(-cpp -fPIC -g)

    # Build type specific flags
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_compile_options(-O0 -check all -traceback -warn all)
    elseif(CMAKE_BUILD_TYPE MATCHES "Release|Fast")
        # Common optimization flags for Release and Fast builds
        add_compile_options(-O3 -fp-model fast=2 -no-prec-div -traceback)

        # Architecture-specific flags
        if(CMAKE_BUILD_TYPE MATCHES Fast)
            add_compile_options(-xHost)  # Native architecture
        else()
            add_compile_options(-xCORE-AVX2)  # Broad compatibility baseline
        endif()
    endif()

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
    # Base flags for all builds
    add_compile_options(-Mpreprocess -fPIC -g)

    # Note: OpenACC is applied selectively to interpolate library only
    # See src/interpolate/CMakeLists.txt

    # Build type specific flags
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_compile_options(-O0 -Mbounds -traceback)
    elseif(CMAKE_BUILD_TYPE MATCHES "Release|Fast")
        # Common optimization flags for Release and Fast builds
        add_compile_options(-fast -O3 -Mfprelaxed -traceback)

        # Architecture-specific flags
        if(CMAKE_BUILD_TYPE MATCHES Fast)
            add_compile_options(-tp=host)  # Native architecture
        endif()
    endif()
endif()

include_directories(${MPI_Fortran_INCLUDE_PATH})
message(STATUS "The MPI_Fortran_INCLUDE_PATH is ${MPI_Fortran_INCLUDE_PATH}")

# Detect and configure external dependencies (HDF5, NetCDF, FFTW)
# Uses smoke tests to verify ABI compatibility with current compiler.
# Falls back to building from source if system libraries are incompatible.
include(${PROJECT_SOURCE_DIR}/cmake/DetectDependencies.cmake)

# Libraries that libneo depends on
add_subdirectory(src/contrib)
add_subdirectory(extra/MyMPILib)
add_subdirectory(src/hdf5_tools)
add_subdirectory(src/polylag)
add_subdirectory(src/interpolate)
add_subdirectory(src/field)
add_subdirectory(src/collisions)
add_subdirectory(src/transport)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Sources for libneo core
include(CMakeSources.in)

# Libraries that depend on libneo
add_subdirectory(src/species)
add_subdirectory(src/magfie)

if(ENABLE_OPENMP)
    target_link_libraries(neo PUBLIC OpenMP::OpenMP_Fortran)
    target_link_libraries(magfie PUBLIC OpenMP::OpenMP_Fortran)
endif()

add_subdirectory(src/efit_to_boozer)
add_executable(efit_to_boozer.x
    app/efit_to_boozer.f90
)
target_link_libraries(efit_to_boozer.x efit_to_boozer)
add_executable(hamada_to_boozer.x
    src/efit_to_boozer/field_line_integration_for_converter.f90
    src/efit_to_boozer/rhs_converter.f90
    app/hamada_to_boozer.f90
)
target_link_libraries(hamada_to_boozer.x efit_to_boozer)

add_executable(vmec_to_chartmap.x
    app/vmec_to_chartmap.f90
)
target_link_libraries(vmec_to_chartmap.x neo)


if(LIBNEO_BUILD_TESTING)
    add_subdirectory(test)
endif()


# MacOS RPATH specifics
# see https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "$ORIGIN/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "$ORIGIN/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "$ORIGIN/lib")
endif("${isSystemDir}" STREQUAL "-1")

# Python bindings

find_package(Python COMPONENTS Interpreter Development NumPy)

message(STATUS "Python: ${Python_EXECUTABLE}")
message(STATUS "Python version: ${Python_VERSION}")
message(STATUS "Python include dirs: ${Python_INCLUDE_DIRS}")
message(STATUS "Python libraries: ${Python_LIBRARIES}")

if(NOT Python_FOUND)
    message("Python with NumPy not found, skipping interface build.")
    return()
endif()

# F2PY headers
execute_process(
    COMMAND "${Python_EXECUTABLE}"
    -c "import os; import numpy.f2py; print(os.path.dirname(numpy.f2py.__file__))"
    OUTPUT_VARIABLE F2PY_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
set(F2PY_INCLUDE_DIR "${F2PY_DIR}/src")

include_directories(
    BEFORE
    ${Python_INCLUDE_DIRS}
    ${Python_NumPy_INCLUDE_DIRS}
    ${F2PY_INCLUDE_DIR}
)

message(STATUS "Python f2py include dir: ${F2PY_INCLUDE_DIR}")
message(STATUS "Python numpy include dir: ${Python_NumPy_INCLUDE_DIRS}")
message(STATUS "Python binary output dir: ${CMAKE_CURRENT_BINARY_DIR}")

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")

add_custom_command(
    OUTPUT _magfiemodule.c _magfie-f2pywrappers.f _magfie-f2pywrappers2.f90
    DEPENDS
        src/f2py_interfaces/f2py_magfie.f90
        src/f2py_interfaces/f2py_vmec_wrappers.f90
        src/f2py_interfaces/f2py_vmec_support.f90
    VERBATIM
    COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
        --f2cmap ${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap
        "${CMAKE_CURRENT_SOURCE_DIR}/src/f2py_interfaces/f2py_magfie.f90"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/f2py_interfaces/f2py_vmec_wrappers.f90"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/f2py_interfaces/f2py_vmec_support.f90"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/magfie/magfie_vmec.f90"
        -m _magfie --lower
)

Python_add_library(_magfie WITH_SOABI
    "${CMAKE_CURRENT_BINARY_DIR}/_magfiemodule.c"
    "${CMAKE_CURRENT_BINARY_DIR}/_magfie-f2pywrappers.f"
    "${CMAKE_CURRENT_BINARY_DIR}/_magfie-f2pywrappers2.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/f2py_interfaces/f2py_magfie.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/f2py_interfaces/f2py_vmec_wrappers.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/f2py_interfaces/f2py_vmec_support.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/magfie/magfie_vmec.f90"
)

target_link_libraries(_magfie PUBLIC fortranobject magfie neo interpolate)
add_dependencies(_magfie neo)
install(TARGETS _magfie DESTINATION .)
