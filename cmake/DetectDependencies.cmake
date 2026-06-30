# DetectDependencies.cmake
# Orchestrates detection of all external dependencies with ABI compatibility checks.
# Automatically falls back to building from source when system libraries are
# incompatible or unavailable.
#
# Usage: include(DetectDependencies) from the main CMakeLists.txt
#
# Options (set before including):
#   PREFER_SYSTEM_LIBS  - Try system libraries first (default: ON)
#   FORCE_FETCH_DEPS    - Always build from source (default: OFF)
#
# After including, the following will be available:
#   - HDF5 targets (hdf5::hdf5_fortran, etc.)
#   - NetCDF linking via NETCDFINCLUDE_DIR and NETCDF_FLIBS_LIST, or netcdf::netcdff

include(${CMAKE_CURRENT_LIST_DIR}/CheckFortranDependency.cmake)

# Check if we must fetch due to compiler incompatibility
should_force_fetch_deps(MUST_FETCH_FORTRAN_DEPS)
if(MUST_FETCH_FORTRAN_DEPS)
    message(STATUS "")
    message(STATUS "=== Dependency Detection ===")
    message(STATUS "NVIDIA/PGI compiler detected - building Fortran dependencies from source")
    message(STATUS "")
endif()

# Track what needs to be built from source
set(NEED_BUILD_HDF5 FALSE)
set(NEED_BUILD_NETCDF FALSE)

#------------------------------------------------------------------------------
# HDF5 Detection
#------------------------------------------------------------------------------
macro(detect_hdf5)
    set(HDF5_USABLE FALSE)

    if(NOT MUST_FETCH_FORTRAN_DEPS AND PREFER_SYSTEM_LIBS)
        find_package(HDF5 QUIET COMPONENTS C HL Fortran Fortran_HL)

        if(HDF5_FOUND)
            message(STATUS "Found system HDF5: ${HDF5_VERSION}")

            check_fortran_dependency(
                "HDF5-Fortran"
                "use hdf5"
                "${HDF5_Fortran_INCLUDE_DIRS}"
                "${HDF5_Fortran_LIBRARIES}"
                HDF5_USABLE
            )

            if(HDF5_USABLE)
                message(STATUS "  -> ABI compatible, using system HDF5")
            else()
                message(STATUS "  -> ABI incompatible, will build from source")
            endif()
        else()
            message(STATUS "HDF5 not found on system")
        endif()
    endif()

    if(NOT HDF5_USABLE)
        set(NEED_BUILD_HDF5 TRUE)
    endif()
endmacro()

#------------------------------------------------------------------------------
# NetCDF Detection
#------------------------------------------------------------------------------
macro(detect_netcdf)
    set(NETCDF_USABLE FALSE)

    if(NOT MUST_FETCH_FORTRAN_DEPS AND PREFER_SYSTEM_LIBS)
        find_program(NF_CONFIG "nf-config")

        if(NF_CONFIG)
            execute_process(COMMAND nf-config --includedir
                OUTPUT_VARIABLE NETCDFINCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
            execute_process(COMMAND nc-config --libdir
                OUTPUT_VARIABLE NETCDFLIB_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
            execute_process(COMMAND nf-config --flibs
                OUTPUT_VARIABLE NETCDF_FLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)

            message(STATUS "Found system NetCDF-Fortran: ${NETCDFINCLUDE_DIR}")

            string(REPLACE " " ";" NETCDF_FLIBS_LIST "${NETCDF_FLIBS}")

            check_fortran_dependency(
                "NetCDF-Fortran"
                "use netcdf"
                "${NETCDFINCLUDE_DIR}"
                "${NETCDF_FLIBS_LIST}"
                NETCDF_USABLE
            )

            if(NETCDF_USABLE)
                message(STATUS "  -> ABI compatible, using system NetCDF")
                set(NETCDFINCLUDE_DIR "${NETCDFINCLUDE_DIR}" CACHE PATH "" FORCE)
                set(NETCDFLIB_DIR "${NETCDFLIB_DIR}" CACHE PATH "" FORCE)
                set(NETCDF_FLIBS_LIST "${NETCDF_FLIBS_LIST}" CACHE STRING "" FORCE)
                set(NETCDF_FOUND TRUE CACHE BOOL "" FORCE)
            else()
                message(STATUS "  -> ABI incompatible, will build from source")
            endif()
        else()
            message(STATUS "NetCDF-Fortran not found on system (nf-config missing)")
        endif()
    endif()

    if(NOT NETCDF_USABLE)
        set(NEED_BUILD_NETCDF TRUE)
        # NetCDF depends on HDF5, so we need to build that too
        set(NEED_BUILD_HDF5 TRUE)
    endif()
endmacro()

#------------------------------------------------------------------------------
# Apply detected NetCDF to build (for system libs case)
#------------------------------------------------------------------------------
macro(apply_netcdf_to_build)
    if(NETCDF_USABLE AND NOT NETCDF_FETCHED)
        include_directories(${NETCDFINCLUDE_DIR})
        link_directories(${NETCDFLIB_DIR})
        add_link_options(${NETCDF_FLIBS_LIST})
    endif()
endmacro()

#------------------------------------------------------------------------------
# Main detection sequence
#------------------------------------------------------------------------------
message(STATUS "")
message(STATUS "=== Detecting External Dependencies ===")
message(STATUS "")

# Light consumers (LIBNEO_BUILD_NUMERICS=OFF) link no HDF5-backed target and get
# NetCDF via a pkg-config imported target (set up in the top-level CMakeLists),
# so skip the nf-config detection here: it never triggers the HDF5/NetCDF
# from-source superbuild and never injects nf-config's full static --flibs
# (-lhdf5 -lz -lcurl, no -L) globally onto consumer executables.
if(NOT DEFINED LIBNEO_BUILD_NUMERICS OR LIBNEO_BUILD_NUMERICS)
    detect_hdf5()
    detect_netcdf()
endif()

# If any dependencies need to be built from source, use the superbuild
if(NEED_BUILD_HDF5 OR NEED_BUILD_NETCDF)
    message(STATUS "")
    message(STATUS "Some dependencies will be built from source:")
    if(NEED_BUILD_HDF5)
        message(STATUS "  - HDF5")
    endif()
    if(NEED_BUILD_NETCDF)
        message(STATUS "  - NetCDF-C and NetCDF-Fortran")
    endif()

    # Include the superbuild module and build all needed dependencies
    include(${CMAKE_CURRENT_LIST_DIR}/BuildExternalDependencies.cmake)
    build_all_external_dependencies()

    # Set include directory for built NetCDF
    if(NEED_BUILD_NETCDF)
        include_directories(${DEPS_PREFIX}/include)
    endif()
else()
    # Apply NetCDF settings to the build (system libs)
    apply_netcdf_to_build()
endif()

message(STATUS "")
message(STATUS "=== Dependency Detection Complete ===")
message(STATUS "")
