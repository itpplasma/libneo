# This file is intended for configuration of libneo.
# List of source files are in CMakeSources.in.

# Basic cmake settings/requirements
cmake_minimum_required(VERSION 3.18)

# Basic project settings.
project(libneo)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# The version number.
set (libneo_VERSION_MAJOR 0)
set (libneo_VERSION_MINOR 0)
set (libneo_VERSION_PATCH 3)

enable_language(Fortran)
enable_testing()

if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/BuildConfig.cmake.in)
  include (${CMAKE_CURRENT_BINARY_DIR}/BuildConfig.cmake.in)
endif()

find_program (BASH_PROGRAM bash)

find_package (MPI REQUIRED COMPONENTS Fortran)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")

  set (CMAKE_Fortran_FLAGS "-cpp -g -Og -fPIC -fcheck=all -fbacktrace \
                            -fno-realloc-lhs -fopenmp -Wall -Wextra")

  # https://github.com/numpy/numpy/issues/25777
  set (CMAKE_C_FLAGS "-Wno-error=incompatible-pointer-types")

elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

  set (CMAKE_Fortran_FLAGS "-cpp -check all -traceback \
                            -assume norealloc_lhs -qopenmp -warn all")

endif ()

include_directories(${MPI_Fortran_INCLUDE_PATH})
message(STATUS "The MPI_Fortran_INCLUDE_PATH is ${MPI_Fortran_INCLUDE_PATH}")

find_program(NF_CONFIG "nf-config")

if (NF_CONFIG)
execute_process(COMMAND nf-config --includedir
                OUTPUT_VARIABLE NETCDFINCLUDE_DIR)
execute_process(COMMAND nc-config --libdir
				OUTPUT_VARIABLE NETCDFLIB_DIR)
execute_process(COMMAND nf-config --flibs
                OUTPUT_VARIABLE NETCDF_FLIBS)
else()
message(SEND_ERROR "nf-config not found. Please install libnetcdff-dev")
endif()

string(STRIP ${NETCDFINCLUDE_DIR} NETCDFINCLUDE_DIR)
string(STRIP ${NETCDFLIB_DIR} NETCDFLIB_DIR)
string(STRIP ${NETCDF_FLIBS} NETCDF_FLIBS)

message(STATUS "NetCDF include path: " ${NETCDFINCLUDE_DIR})
message(STATUS "NetCDF lib path: " ${NETCDFLIB_DIR})
message(STATUS "NetCDF Fortran libs: " ${NETCDF_FLIBS})

# Replace space by semicolon in the Fortran libs
string(REPLACE " " ";" NETCDF_FLIBS ${NETCDF_FLIBS})

include_directories(${NETCDFINCLUDE_DIR})
link_directories(${NETCDFLIB_DIR})
add_link_options(${NETCDF_FLIBS})

add_subdirectory(src/MyMPILib)
add_subdirectory(src/hdf5_tools)
add_subdirectory(src/polylag)
add_subdirectory(src/interpolate)
add_subdirectory(src/field)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

include(CMakeSources.in)

add_subdirectory(test)


# MacOS RPATH specifics
# see https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

# Python bindings

# Grab Python
find_package(
  Python
  COMPONENTS Interpreter Development.Module numpy
)

if (Python_FOUND)
    message(STATUS "Python found: ${Python_EXECUTABLE}")
    message(STATUS "Python version: ${Python_VERSION}")
    message(STATUS "Python include dirs: ${Python_INCLUDE_DIRS}")
    message(STATUS "Python libraries: ${Python_LIBRARIES}")
else()
    message("Python with NumPy not found, skipping interface build.")
    return()
endif()

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

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

add_custom_command(
  OUTPUT _magfiemodule.c _magfie-f2pywrappers.f
  DEPENDS src/f2py_interfaces/f2py_magfie.f90
  VERBATIM
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
  "${CMAKE_CURRENT_SOURCE_DIR}/src/f2py_interfaces/f2py_magfie.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/magfie/magfie_vmec.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/spline_vmec_data.f90"
  -m _magfie --lower
)

Python_add_library(_magfie WITH_SOABI
  "${CMAKE_CURRENT_BINARY_DIR}/_magfiemodule.c"
  "${CMAKE_CURRENT_BINARY_DIR}/_magfie-f2pywrappers.f"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/f2py_interfaces/f2py_magfie.f90"
  )

target_link_libraries(_magfie PUBLIC fortranobject magfie)
install(TARGETS _magfie DESTINATION .)

# Check if f90wrap package ist installed
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m "f90wrap"
    RESULT_VARIABLE F90WRAP_RESULT
    OUTPUT_QUIET
    ERROR_QUIET
)

if (F90WRAP_RESULT EQUAL 0)
    message(STATUS "Python f90wrap found.")
else()
    message(STATUS "Python f90wrap not found, skipping further interface builds.")
    return()
endif()

set(FILES_TO_WRAP
  interpolate/interpolate.f90
)
list(TRANSFORM FILES_TO_WRAP PREPEND ${PROJECT_SOURCE_DIR}/src/)

foreach(file ${FILES_TO_WRAP})
    get_filename_component(basename ${file} NAME_WE)
    list(APPEND PREPROCESSED_SOURCES
        ${CMAKE_CURRENT_BINARY_DIR}/${basename}.f90.i
    )
    list(APPEND F90WRAP_OUTPUT
        ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_${basename}.f90.f90
    )
endforeach()

if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_toplevel.f90)
  list(APPEND F90WRAP_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_toplevel.f90)
endif()

list(REMOVE_ITEM FILES_TO_WRAP ${PROJECT_SOURCE_DIR}/src/interpolate.f90)
list(APPEND FILES_TO_WRAP ${PROJECT_SOURCE_DIR}/src/interpolate/interpolate.f90)

set(F2PY_MODULE_NAME "_interpolate")

set(F2PY_MODULE_C "${F2PY_MODULE_NAME}module.c")
set(GENERATED_MODULE_FILE "${F2PY_MODULE_NAME}.${Python_SOABI}")
set(F2PY_WRAPPERS "${CMAKE_CURRENT_BINARY_DIR}/${F2PY_MODULE_C}")

add_custom_target(F90WrapOutput)

foreach (file_to_wrap ${FILES_TO_WRAP})
    get_filename_component(basename ${file_to_wrap} NAME_WE)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${basename}.f90.i
        COMMAND gfortran -E -cpp
            -o ${CMAKE_CURRENT_BINARY_DIR}/${basename}.f90.i
            ${file_to_wrap}
    DEPENDS ${file_to_wrap}
)
endforeach()

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/interpolate.py ${F90WRAP_OUTPUT}
	COMMAND ${Python_EXECUTABLE}
        -m f90wrap
        --f90wrap
        -m interpolate
        -k ${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap
        ${CMAKE_CURRENT_BINARY_DIR}/*.f90.i
        > ${CMAKE_CURRENT_BINARY_DIR}/f90wrap.log 2>&1
	DEPENDS ${PREPROCESSED_SOURCES}
    COMMENT "Processing preprocessed sources with f90wrap"
)

add_custom_command(OUTPUT ${F2PY_WRAPPERS}
    COMMAND ${Python_EXECUTABLE}
        -m f90wrap
        --f2py-f90wrap
        ${F90WRAP_OUTPUT}
        -m ${F2PY_MODULE_NAME}
        --f2cmap ${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap
        --lower
        > ${CMAKE_CURRENT_BINARY_DIR}/f2py.log 2>&1
    DEPENDS ${F90WRAP_OUTPUT}
    COMMENT "Processing wrapped sources with f2py"
)

Python_add_library(${F2PY_MODULE_NAME} WITH_SOABI
   ${F2PY_WRAPPERS}
   ${F2PY_INCLUDE_DIR}/fortranobject.c
   ${F90WRAP_OUTPUT}
)

add_dependencies(${F2PY_MODULE_NAME} neo)

target_link_libraries(${F2PY_MODULE_NAME} PRIVATE neo)

set_target_properties(
    ${F2PY_MODULE_NAME}
    PROPERTIES
    PREFIX ""
    OUTPUT_NAME "${GENERATED_MODULE_FILE}"
    SUFFIX ".so"
    LINKER_LANGUAGE C
)

install(TARGETS ${F2PY_MODULE_NAME} DESTINATION .)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/interpolate.py DESTINATION .)
