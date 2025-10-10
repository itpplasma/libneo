# This file is intended for configuration of libneo.
# List of source files are in CMakeSources.in.

# Basic cmake settings/requirements
cmake_minimum_required(VERSION 3.22)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Basic project settings.
project(efit_to_boozer)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# The version number.
set (efit_to_boozer_VERSION_MAJOR 0)
set (efit_to_boozer_VERSION_MINOR 0)
set (efit_to_boozer_VERSION_PATCH 3)

enable_language(C Fortran)

if(LIBNEO_ENABLE_TESTS)
    enable_testing()
endif()

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

  set (CMAKE_Fortran_FLAGS "-cpp -g -Og -fcheck=all -fbacktrace \
                            -fno-realloc-lhs -Wall -Wextra")

  # https://github.com/numpy/numpy/issues/25777
  set (CMAKE_C_FLAGS "-Wno-error=incompatible-pointer-types")

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

  set (CMAKE_Fortran_FLAGS "-cpp -check all -traceback \
                            -assume norealloc_lhs -qopenmp -warn all")

endif ()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# MacOS RPATH specifics
# see https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "$ORIGIN")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "$ORIGIN" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "$ORIGIN")
endif("${isSystemDir}" STREQUAL "-1")

add_library(efit_to_boozer STATIC
    efit_to_boozer_mod.f90
    field_line_integration_for_Boozer.f90
    rhs.f90
    spline_and_interpolate_magdata.f90
)
target_link_libraries(efit_to_boozer
    BLAS::BLAS
    LAPACK::LAPACK
    neo
    magfie
)

# Grab Python
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

message(STATUS "Python found: ${Python_EXECUTABLE}")
message(STATUS "Python version: ${Python_VERSION}")
message(STATUS "Python include dirs: ${Python_INCLUDE_DIRS}")
message(STATUS "Python libraries: ${Python_LIBRARIES}")

# F2PY headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# Grab the variables from a local Python installation
# F2PY headers
execute_process(
    COMMAND "${Python_EXECUTABLE}"
    -c "import numpy; print(numpy.__version__)"
    OUTPUT_VARIABLE NUMPY_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if (${NUMPY_VERSION} VERSION_GREATER_EQUAL "1.21.1")
	execute_process(
        COMMAND "${Python_EXECUTABLE}"
        -c "import numpy.f2py; print(numpy.f2py.get_include())"
        OUTPUT_VARIABLE F2PY_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
	)
else()
    execute_process(
        COMMAND "${Python_EXECUTABLE}"
        -c "import numpy; print(numpy.get_include())"
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(F2PY_INCLUDE_DIR "${NUMPY_INCLUDE_DIR}/../../f2py/src/")
endif()

include_directories(
    BEFORE
    ${Python_INCLUDE_DIRS}
    ${Python_NumPy_INCLUDE_DIRS}
    ${F2PY_INCLUDE_DIR}
)

message(STATUS "Python f2py include dir: ${F2PY_INCLUDE_DIR}")
message(STATUS "Python numpy include dir: ${Python_NumPy_INCLUDE_DIRS}")
message(STATUS "Python binary output dir: ${CMAKE_CURRENT_BINARY_DIR}")

add_custom_command(
  OUTPUT _efit_to_boozermodule.c _efit_to_boozer-f2pywrappers.f _efit_to_boozer-f2pywrappers2.f90
  DEPENDS f2py_interfaces/f2py_efit_to_boozer.f90
  VERBATIM
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
  "${CMAKE_CURRENT_SOURCE_DIR}/f2py_interfaces/f2py_efit_to_boozer.f90"
  "${CMAKE_SOURCE_DIR}/src/magfie/field_divB0.f90"
  "${CMAKE_SOURCE_DIR}/src/magfie/field_eq_mod.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/efit_to_boozer_mod.f90"
  "${CMAKE_SOURCE_DIR}/src/spl_three_to_five.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/field_line_integration_for_Boozer.f90"
  "${CMAKE_SOURCE_DIR}/src/plag_coeff.f90"
  "${CMAKE_SOURCE_DIR}/src/binsrc.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/rhs.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/spline_and_interpolate_magdata.f90"
  -m _efit_to_boozer --lower skip: oddorderspline polleg binomial odeint_allroutines', alloc_odeint odeint rkck rkqs :
)

Python_add_library(_efit_to_boozer WITH_SOABI
  "${CMAKE_CURRENT_BINARY_DIR}/_efit_to_boozermodule.c"
  "${CMAKE_CURRENT_BINARY_DIR}/_efit_to_boozer-f2pywrappers.f"
  "${CMAKE_CURRENT_BINARY_DIR}/_efit_to_boozer-f2pywrappers2.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/f2py_interfaces/f2py_efit_to_boozer.f90"
)

target_link_libraries(_efit_to_boozer PUBLIC fortranobject efit_to_boozer)
install(TARGETS _efit_to_boozer DESTINATION .)
