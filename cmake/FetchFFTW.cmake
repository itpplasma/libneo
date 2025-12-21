# FetchFFTW.cmake
# Downloads and builds FFTW3 from source.

include(FetchContent)

message(STATUS "Building FFTW from source...")

# FFTW 3.3.10 - latest stable
FetchContent_Declare(fftw3
    URL https://www.fftw.org/fftw-3.3.10.tar.gz
    URL_HASH SHA256=56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# FFTW uses autotools, but also has CMake support
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(ENABLE_OPENMP ON CACHE BOOL "" FORCE)
set(ENABLE_THREADS ON CACHE BOOL "" FORCE)
set(ENABLE_FLOAT OFF CACHE BOOL "" FORCE)
set(ENABLE_LONG_DOUBLE OFF CACHE BOOL "" FORCE)
set(ENABLE_QUAD_PRECISION OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(fftw3)

# Create aliases to match FindFFTW target names
if(TARGET fftw3)
    add_library(FFTW::Double ALIAS fftw3)
endif()
if(TARGET fftw3_omp)
    add_library(FFTW::DoubleOpenMP ALIAS fftw3_omp)
elseif(TARGET fftw3_threads)
    add_library(FFTW::DoubleThreads ALIAS fftw3_threads)
endif()

set(FFTW_FOUND TRUE CACHE BOOL "" FORCE)
set(FFTW_FETCHED TRUE CACHE BOOL "" FORCE)

message(STATUS "FFTW will be built from source")
