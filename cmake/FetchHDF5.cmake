# FetchHDF5.cmake
# Downloads and builds HDF5 with Fortran support from source.

include(FetchContent)

message(STATUS "Building HDF5 from source...")

# HDF5 1.14.5 - latest stable with good Fortran support
FetchContent_Declare(hdf5
    URL https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.5/hdf5-1.14.5.tar.gz
    URL_HASH SHA256=ec2e13c52571f997d4fa98893d59fef2cbc1a53c127712e09e0c3f3e4e4b2bbb
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# Configure HDF5 build options
set(HDF5_BUILD_FORTRAN ON CACHE BOOL "" FORCE)
set(HDF5_BUILD_CPP_LIB OFF CACHE BOOL "" FORCE)
set(HDF5_BUILD_JAVA OFF CACHE BOOL "" FORCE)
set(HDF5_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(HDF5_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
set(HDF5_BUILD_HL_LIB ON CACHE BOOL "" FORCE)
set(HDF5_ENABLE_PARALLEL OFF CACHE BOOL "" FORCE)
set(HDF5_ENABLE_Z_LIB_SUPPORT ON CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(HDF5_BUILD_TESTING OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(hdf5)

# Create aliases to match find_package(HDF5) target names
if(TARGET hdf5-static)
    add_library(hdf5::hdf5 ALIAS hdf5-static)
endif()
if(TARGET hdf5_hl-static)
    add_library(hdf5::hdf5_hl ALIAS hdf5_hl-static)
endif()
if(TARGET hdf5_fortran-static)
    add_library(hdf5::hdf5_fortran ALIAS hdf5_fortran-static)
endif()
if(TARGET hdf5_hl_fortran-static)
    add_library(hdf5::hdf5_hl_fortran ALIAS hdf5_hl_fortran-static)
endif()

set(HDF5_FOUND TRUE CACHE BOOL "" FORCE)
set(HDF5_FETCHED TRUE CACHE BOOL "" FORCE)

message(STATUS "HDF5 will be built from source")
