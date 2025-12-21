# FetchNetCDF.cmake
# Downloads and builds NetCDF-C and NetCDF-Fortran from source.
# Requires HDF5 to be available (either found or fetched).

include(FetchContent)

message(STATUS "Building NetCDF from source...")

# NetCDF-C 4.9.2
FetchContent_Declare(netcdf-c
    GIT_REPOSITORY https://github.com/Unidata/netcdf-c.git
    GIT_TAG v4.9.2
    GIT_SHALLOW TRUE
)

# Configure NetCDF-C build options - minimal build
set(ENABLE_DAP OFF CACHE BOOL "" FORCE)
set(ENABLE_BYTERANGE OFF CACHE BOOL "" FORCE)
set(ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(ENABLE_EXTREME_NUMBERS OFF CACHE BOOL "" FORCE)
set(ENABLE_PARALLEL4 OFF CACHE BOOL "" FORCE)
set(ENABLE_PNETCDF OFF CACHE BOOL "" FORCE)
set(ENABLE_CDF5 ON CACHE BOOL "" FORCE)
set(ENABLE_NCZARR OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_UTILITIES OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(netcdf-c)

# Get NetCDF-C paths for NetCDF-Fortran
FetchContent_GetProperties(netcdf-c)
set(NETCDF_C_INCLUDE_DIR "${netcdf-c_BINARY_DIR}/include" CACHE PATH "" FORCE)
set(NETCDF_C_LIBRARY "${netcdf-c_BINARY_DIR}/liblib/libnetcdf.a" CACHE FILEPATH "" FORCE)

# NetCDF-Fortran 4.6.1
FetchContent_Declare(netcdf-fortran
    GIT_REPOSITORY https://github.com/Unidata/netcdf-fortran.git
    GIT_TAG v4.6.1
    GIT_SHALLOW TRUE
)

# Configure NetCDF-Fortran build options
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(netcdf-fortran)

# Create targets that match what nf-config would provide
FetchContent_GetProperties(netcdf-fortran)
set(NETCDF_FORTRAN_INCLUDE_DIR "${netcdf-fortran_BINARY_DIR}/fortran" CACHE PATH "" FORCE)

# Set variables for downstream use
set(NETCDF_FOUND TRUE CACHE BOOL "" FORCE)
set(NETCDF_FETCHED TRUE CACHE BOOL "" FORCE)
set(NETCDFINCLUDE_DIR "${netcdf-fortran_BINARY_DIR}/fortran" CACHE PATH "" FORCE)

message(STATUS "NetCDF-C and NetCDF-Fortran will be built from source")
