# BuildExternalDependencies.cmake
# Superbuild pattern for external Fortran dependencies.
#
# Uses ExternalProject_Add to build dependencies at BUILD time (not configure time).
# This ensures proper sequencing: HDF5 is installed before NetCDF-C configures,
# NetCDF-C is installed before NetCDF-Fortran configures, etc.
#
# The dependencies are installed to a common prefix that the main build can find.

include(ExternalProject)

set(DEPS_PREFIX "${CMAKE_BINARY_DIR}/deps" CACHE PATH "Installation prefix for dependencies")
set(DEPS_SOURCE_DIR "${CMAKE_BINARY_DIR}/deps-src" CACHE PATH "Source directory for dependencies")
set(DEPS_BUILD_DIR "${CMAKE_BINARY_DIR}/deps-build" CACHE PATH "Build directory for dependencies")

# Create the prefix directory structure
file(MAKE_DIRECTORY ${DEPS_PREFIX}/lib)
file(MAKE_DIRECTORY ${DEPS_PREFIX}/include)

# Common CMake arguments for all dependencies
# -fPIC is needed so static libs can be linked into shared objects (Python modules)
set(COMMON_CMAKE_ARGS
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${DEPS_PREFIX}
    -DCMAKE_PREFIX_PATH=${DEPS_PREFIX}
    -DBUILD_SHARED_LIBS=OFF
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    -DCMAKE_C_FLAGS=-fPIC
    -DCMAKE_Fortran_FLAGS=-fPIC
)

#------------------------------------------------------------------------------
# HDF5
#------------------------------------------------------------------------------
function(build_hdf5)
    message(STATUS "Will build HDF5 1.14.5 from source")

    ExternalProject_Add(hdf5_external
        URL https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.5/hdf5-1.14.5.tar.gz
        URL_HASH SHA256=ec2e13c52e60f9a01491bb3158cb3778c985697131fc6a342262d32a26e58e44
        DOWNLOAD_DIR ${DEPS_SOURCE_DIR}
        SOURCE_DIR ${DEPS_SOURCE_DIR}/hdf5
        BINARY_DIR ${DEPS_BUILD_DIR}/hdf5
        INSTALL_DIR ${DEPS_PREFIX}
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        CMAKE_ARGS
            ${COMMON_CMAKE_ARGS}
            -DHDF5_BUILD_FORTRAN=ON
            -DHDF5_BUILD_CPP_LIB=OFF
            -DHDF5_BUILD_JAVA=OFF
            -DHDF5_BUILD_EXAMPLES=OFF
            -DHDF5_BUILD_TOOLS=OFF
            -DHDF5_BUILD_HL_LIB=ON
            -DHDF5_ENABLE_PARALLEL=OFF
            -DHDF5_ENABLE_Z_LIB_SUPPORT=ON
            -DBUILD_TESTING=OFF
            -DHDF5_BUILD_TESTING=OFF
        BUILD_BYPRODUCTS
            ${DEPS_PREFIX}/lib/libhdf5.a
            ${DEPS_PREFIX}/lib/libhdf5_hl.a
            ${DEPS_PREFIX}/lib/libhdf5_f90cstub.a
            ${DEPS_PREFIX}/lib/libhdf5_hl_f90cstub.a
            ${DEPS_PREFIX}/lib/libhdf5_fortran.a
            ${DEPS_PREFIX}/lib/libhdf5_hl_fortran.a
    )

    # Create imported targets that the main build can use
    # HDF5 creates separate C stub libraries for Fortran bindings
    # Use if(NOT TARGET) guards in case find_package already created these targets
    if(NOT TARGET hdf5::hdf5)
        add_library(hdf5::hdf5 STATIC IMPORTED GLOBAL)
    endif()
    if(NOT TARGET hdf5::hdf5_hl)
        add_library(hdf5::hdf5_hl STATIC IMPORTED GLOBAL)
    endif()
    if(NOT TARGET hdf5::hdf5_f90cstub)
        add_library(hdf5::hdf5_f90cstub STATIC IMPORTED GLOBAL)
    endif()
    if(NOT TARGET hdf5::hdf5_hl_f90cstub)
        add_library(hdf5::hdf5_hl_f90cstub STATIC IMPORTED GLOBAL)
    endif()
    if(NOT TARGET hdf5::hdf5_fortran)
        add_library(hdf5::hdf5_fortran STATIC IMPORTED GLOBAL)
    endif()
    if(NOT TARGET hdf5::hdf5_hl_fortran)
        add_library(hdf5::hdf5_hl_fortran STATIC IMPORTED GLOBAL)
    endif()

    # For static linking, specify all dependencies explicitly to ensure correct link order
    set_target_properties(hdf5::hdf5 PROPERTIES
        IMPORTED_LOCATION ${DEPS_PREFIX}/lib/libhdf5.a
        INTERFACE_INCLUDE_DIRECTORIES ${DEPS_PREFIX}/include
        INTERFACE_LINK_LIBRARIES "z;dl;m"
    )
    set_target_properties(hdf5::hdf5_hl PROPERTIES
        IMPORTED_LOCATION ${DEPS_PREFIX}/lib/libhdf5_hl.a
        INTERFACE_INCLUDE_DIRECTORIES ${DEPS_PREFIX}/include
        INTERFACE_LINK_LIBRARIES "hdf5::hdf5;z;dl;m"
    )
    set_target_properties(hdf5::hdf5_f90cstub PROPERTIES
        IMPORTED_LOCATION ${DEPS_PREFIX}/lib/libhdf5_f90cstub.a
        INTERFACE_INCLUDE_DIRECTORIES ${DEPS_PREFIX}/include
        INTERFACE_LINK_LIBRARIES "hdf5::hdf5;z;dl;m"
    )
    set_target_properties(hdf5::hdf5_hl_f90cstub PROPERTIES
        IMPORTED_LOCATION ${DEPS_PREFIX}/lib/libhdf5_hl_f90cstub.a
        INTERFACE_INCLUDE_DIRECTORIES ${DEPS_PREFIX}/include
        INTERFACE_LINK_LIBRARIES "hdf5::hdf5_hl;hdf5::hdf5;z;dl;m"
    )
    set_target_properties(hdf5::hdf5_fortran PROPERTIES
        IMPORTED_LOCATION ${DEPS_PREFIX}/lib/libhdf5_fortran.a
        INTERFACE_INCLUDE_DIRECTORIES ${DEPS_PREFIX}/include
        INTERFACE_LINK_LIBRARIES "hdf5::hdf5_f90cstub;hdf5::hdf5;z;dl;m"
    )
    set_target_properties(hdf5::hdf5_hl_fortran PROPERTIES
        IMPORTED_LOCATION ${DEPS_PREFIX}/lib/libhdf5_hl_fortran.a
        INTERFACE_INCLUDE_DIRECTORIES ${DEPS_PREFIX}/include
        INTERFACE_LINK_LIBRARIES "hdf5::hdf5_hl_f90cstub;hdf5::hdf5_fortran;hdf5::hdf5_f90cstub;hdf5::hdf5_hl;hdf5::hdf5;z;dl;m"
    )

    # These targets depend on hdf5_external being built
    add_dependencies(hdf5::hdf5 hdf5_external)
    add_dependencies(hdf5::hdf5_hl hdf5_external)
    add_dependencies(hdf5::hdf5_f90cstub hdf5_external)
    add_dependencies(hdf5::hdf5_hl_f90cstub hdf5_external)
    add_dependencies(hdf5::hdf5_fortran hdf5_external)
    add_dependencies(hdf5::hdf5_hl_fortran hdf5_external)

    set(HDF5_FOUND TRUE CACHE BOOL "" FORCE)
    set(HDF5_FETCHED TRUE CACHE BOOL "" FORCE)
endfunction()

#------------------------------------------------------------------------------
# NetCDF-C
#------------------------------------------------------------------------------
function(build_netcdf_c)
    message(STATUS "Will build NetCDF-C 4.9.2 from source")

    ExternalProject_Add(netcdf_c_external
        GIT_REPOSITORY https://github.com/Unidata/netcdf-c.git
        GIT_TAG v4.9.2
        GIT_SHALLOW TRUE
        DOWNLOAD_DIR ${DEPS_SOURCE_DIR}
        SOURCE_DIR ${DEPS_SOURCE_DIR}/netcdf-c
        BINARY_DIR ${DEPS_BUILD_DIR}/netcdf-c
        INSTALL_DIR ${DEPS_PREFIX}
        CMAKE_ARGS
            ${COMMON_CMAKE_ARGS}
            -DCMAKE_PREFIX_PATH=${DEPS_PREFIX}
            -DHDF5_ROOT=${DEPS_PREFIX}
            -DHDF5_DIR=${DEPS_PREFIX}
            -DENABLE_DAP=OFF
            -DENABLE_BYTERANGE=OFF
            -DENABLE_EXAMPLES=OFF
            -DENABLE_TESTS=OFF
            -DENABLE_EXTREME_NUMBERS=OFF
            -DENABLE_PARALLEL4=OFF
            -DENABLE_PNETCDF=OFF
            -DENABLE_CDF5=ON
            -DENABLE_NCZARR=OFF
            -DBUILD_UTILITIES=OFF
            -DENABLE_FILTER_TESTING=OFF
            -DENABLE_PLUGINS=OFF
            -DENABLE_HDF5=ON
            -DENABLE_NETCDF_4=ON
        DEPENDS hdf5_external
        BUILD_BYPRODUCTS ${DEPS_PREFIX}/lib/libnetcdf.a
    )

    add_library(netcdf::netcdf STATIC IMPORTED GLOBAL)
    set_target_properties(netcdf::netcdf PROPERTIES
        IMPORTED_LOCATION ${DEPS_PREFIX}/lib/libnetcdf.a
        INTERFACE_INCLUDE_DIRECTORIES ${DEPS_PREFIX}/include
        INTERFACE_LINK_LIBRARIES "hdf5::hdf5_hl;hdf5::hdf5;z"
    )
    add_dependencies(netcdf::netcdf netcdf_c_external)
endfunction()

#------------------------------------------------------------------------------
# NetCDF-Fortran
#------------------------------------------------------------------------------
function(build_netcdf_fortran)
    message(STATUS "Will build NetCDF-Fortran 4.6.1 from source")

    ExternalProject_Add(netcdf_fortran_external
        GIT_REPOSITORY https://github.com/Unidata/netcdf-fortran.git
        GIT_TAG v4.6.1
        GIT_SHALLOW TRUE
        DOWNLOAD_DIR ${DEPS_SOURCE_DIR}
        SOURCE_DIR ${DEPS_SOURCE_DIR}/netcdf-fortran
        BINARY_DIR ${DEPS_BUILD_DIR}/netcdf-fortran
        INSTALL_DIR ${DEPS_PREFIX}
        CMAKE_ARGS
            ${COMMON_CMAKE_ARGS}
            -DCMAKE_PREFIX_PATH=${DEPS_PREFIX}
            -DNETCDF_C_LIBRARY=${DEPS_PREFIX}/lib/libnetcdf.a
            -DNETCDF_C_INCLUDE_DIR=${DEPS_PREFIX}/include
            -DENABLE_TESTS=OFF
            -DENABLE_EXAMPLES=OFF
        DEPENDS netcdf_c_external
        BUILD_BYPRODUCTS ${DEPS_PREFIX}/lib/libnetcdff.a
    )

    add_library(netcdf::netcdff STATIC IMPORTED GLOBAL)
    set_target_properties(netcdf::netcdff PROPERTIES
        IMPORTED_LOCATION ${DEPS_PREFIX}/lib/libnetcdff.a
        INTERFACE_INCLUDE_DIRECTORIES ${DEPS_PREFIX}/include
        INTERFACE_LINK_LIBRARIES netcdf::netcdf
    )
    add_dependencies(netcdf::netcdff netcdf_fortran_external)

    # Set variables for legacy code that uses NETCDFINCLUDE_DIR
    set(NETCDFINCLUDE_DIR "${DEPS_PREFIX}/include" CACHE PATH "" FORCE)
    set(NETCDF_FLIBS_LIST "netcdf::netcdff" CACHE STRING "" FORCE)
    set(NETCDF_FOUND TRUE CACHE BOOL "" FORCE)
    set(NETCDF_FETCHED TRUE CACHE BOOL "" FORCE)
endfunction()

#------------------------------------------------------------------------------
# FFTW (C library - no Fortran ABI issues, but include for completeness)
#------------------------------------------------------------------------------
function(build_fftw)
    message(STATUS "Will build FFTW 3.3.10 from source")

    ExternalProject_Add(fftw_external
        URL http://www.fftw.org/fftw-3.3.10.tar.gz
        URL_HASH MD5=8ccbf6a5ea78a16dbc3e1306e234cc5c
        DOWNLOAD_DIR ${DEPS_SOURCE_DIR}
        SOURCE_DIR ${DEPS_SOURCE_DIR}/fftw
        BINARY_DIR ${DEPS_BUILD_DIR}/fftw
        INSTALL_DIR ${DEPS_PREFIX}
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        CMAKE_ARGS
            ${COMMON_CMAKE_ARGS}
            -DCMAKE_POLICY_VERSION_MINIMUM=3.5
            -DENABLE_THREADS=ON
            -DENABLE_OPENMP=ON
            -DBUILD_TESTS=OFF
        BUILD_BYPRODUCTS
            ${DEPS_PREFIX}/lib/libfftw3.a
            ${DEPS_PREFIX}/lib/libfftw3_threads.a
    )

    add_library(FFTW::Double STATIC IMPORTED GLOBAL)
    add_library(FFTW::DoubleThreads STATIC IMPORTED GLOBAL)

    set_target_properties(FFTW::Double PROPERTIES
        IMPORTED_LOCATION ${DEPS_PREFIX}/lib/libfftw3.a
        INTERFACE_INCLUDE_DIRECTORIES ${DEPS_PREFIX}/include
    )
    set_target_properties(FFTW::DoubleThreads PROPERTIES
        IMPORTED_LOCATION ${DEPS_PREFIX}/lib/libfftw3_threads.a
        INTERFACE_INCLUDE_DIRECTORIES ${DEPS_PREFIX}/include
        INTERFACE_LINK_LIBRARIES FFTW::Double
    )

    add_dependencies(FFTW::Double fftw_external)
    add_dependencies(FFTW::DoubleThreads fftw_external)

    set(FFTW_FOUND TRUE CACHE BOOL "" FORCE)
    set(FFTW_FETCHED TRUE CACHE BOOL "" FORCE)
endfunction()

#------------------------------------------------------------------------------
# Build all dependencies in order (only those that are needed)
#------------------------------------------------------------------------------
function(build_all_external_dependencies)
    message(STATUS "")
    message(STATUS "=== Building External Dependencies from Source ===")
    message(STATUS "Installation prefix: ${DEPS_PREFIX}")
    message(STATUS "")

    if(NEED_BUILD_HDF5)
        build_hdf5()
    endif()
    if(NEED_BUILD_NETCDF)
        build_netcdf_c()
        build_netcdf_fortran()
    endif()
    if(NEED_BUILD_FFTW)
        build_fftw()
    endif()

    message(STATUS "")
    message(STATUS "=== External Dependencies Configured ===")
    message(STATUS "")
endfunction()
