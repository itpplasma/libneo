# SetupNVHPCCuda.cmake
# Auto-detect CUDA toolkit for NVHPC compilers
#
# NVHPC compilers require NVHPC_CUDA_HOME to be set when using OpenACC GPU offloading.
# This module auto-detects common CUDA installation paths and sets the environment.
#
# Search order:
# 1. User-provided NVHPC_CUDA_HOME (respected if set)
# 2. /opt/cuda (common system installation)
# 3. /usr/local/cuda
# 4. CUDA_HOME environment variable
#
# For MPI wrappers, NVCOMPILER_COMM_LIBS_HOME may also be needed.

function(setup_nvhpc_cuda)
    # Check if we're using NVHPC compiler
    # This runs before compiler detection, so we check the compiler path
    if(CMAKE_Fortran_COMPILER)
        get_filename_component(_fc_name "${CMAKE_Fortran_COMPILER}" NAME)
        string(FIND "${_fc_name}" "nvfortran" _is_nvfortran)
        string(FIND "${_fc_name}" "mpifort" _is_mpifort)
        string(FIND "${CMAKE_Fortran_COMPILER}" "hpc_sdk" _is_hpcsdk)

        # Only proceed if this looks like NVHPC
        if(_is_nvfortran LESS 0 AND NOT (_is_mpifort GREATER_EQUAL 0 AND _is_hpcsdk GREATER_EQUAL 0))
            return()
        endif()
    else()
        return()
    endif()

    # Check if NVHPC_CUDA_HOME is already set
    if(DEFINED ENV{NVHPC_CUDA_HOME})
        message(STATUS "NVHPC_CUDA_HOME already set: $ENV{NVHPC_CUDA_HOME}")
        return()
    endif()

    # Search for CUDA installation
    set(_cuda_search_paths
        "/opt/cuda"
        "/usr/local/cuda"
        "$ENV{CUDA_HOME}"
    )

    set(_cuda_found OFF)
    foreach(_cuda_path ${_cuda_search_paths})
        if(_cuda_path AND EXISTS "${_cuda_path}/bin/nvcc")
            set(ENV{NVHPC_CUDA_HOME} "${_cuda_path}")
            message(STATUS "Auto-detected CUDA for NVHPC: ${_cuda_path}")
            set(_cuda_found ON)
            break()
        endif()
    endforeach()

    if(NOT _cuda_found)
        message(WARNING
            "NVHPC compiler detected but no CUDA toolkit found.\n"
            "GPU offloading may fail. Set NVHPC_CUDA_HOME manually or install CUDA to:\n"
            "  /opt/cuda\n"
            "  /usr/local/cuda")
    endif()

    # For NVHPC MPI wrappers, also set NVCOMPILER_COMM_LIBS_HOME if needed
    if(_is_hpcsdk GREATER_EQUAL 0 AND _is_mpifort GREATER_EQUAL 0)
        if(NOT DEFINED ENV{NVCOMPILER_COMM_LIBS_HOME})
            # Extract the HPC SDK path from the compiler path
            string(REGEX REPLACE "/comm_libs/.*" "/comm_libs" _comm_libs_path "${CMAKE_Fortran_COMPILER}")
            if(EXISTS "${_comm_libs_path}")
                set(ENV{NVCOMPILER_COMM_LIBS_HOME} "${_comm_libs_path}")
                message(STATUS "Auto-detected NVCOMPILER_COMM_LIBS_HOME: ${_comm_libs_path}")
            endif()
        endif()
    endif()
endfunction()

# Call the function when this module is included
setup_nvhpc_cuda()
