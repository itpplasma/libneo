cmake_minimum_required(VERSION 3.22)

project(libneo_gpu_spline_1d_many LANGUAGES C Fortran)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_C_STANDARD 11)

option(DRAFT_ENABLE_OPENACC "Build OpenACC variant (NVHPC)" ON)
option(DRAFT_ENABLE_OMP_TARGET "Build OpenMP target offload variant (NVHPC)" ON)
option(DRAFT_ENABLE_CUDA_FORTRAN "Build CUDA Fortran variant (NVHPC)" ON)
option(DRAFT_ENABLE_CUDA_C "Build CUDA C kernel variant (NVCC + Fortran)" ON)

if(NOT CMAKE_Fortran_COMPILER MATCHES "nvfortran")
    message(FATAL_ERROR "This draft benchmark expects NVHPC nvfortran as Fortran compiler.")
endif()

set(DRAFT_ROOT ${CMAKE_CURRENT_LIST_DIR})
set(LIBNEO_ROOT ${DRAFT_ROOT}/../..)

add_library(draft_interpolate STATIC
    ${LIBNEO_ROOT}/src/spl_three_to_five.f90
    ${LIBNEO_ROOT}/src/interpolate/batch_interpolate_types.f90
    ${LIBNEO_ROOT}/src/interpolate/batch_interpolate_1d.f90
)
target_include_directories(draft_interpolate PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

add_library(draft_spline1d_many STATIC
    src/spline1d_many_cpu.f90
    src/spline1d_many_openacc.f90
    src/spline1d_many_omp_target.f90
    src/spline1d_many_cuda_fortran.f90
    src/spline1d_many_cuda_c_wrapper.f90
    src/util_lcg_rng.f90
)
target_link_libraries(draft_spline1d_many PUBLIC draft_interpolate)

add_executable(bench_spline1d_many src/bench_spline1d_many.f90)
target_link_libraries(bench_spline1d_many PRIVATE draft_spline1d_many)

target_link_options(bench_spline1d_many PRIVATE
    -gpu=ccnative
)

if(DRAFT_ENABLE_OPENACC)
    target_link_options(bench_spline1d_many PRIVATE -acc)
endif()

if(DRAFT_ENABLE_OMP_TARGET)
    target_link_options(bench_spline1d_many PRIVATE -mp=gpu)
endif()

if(DRAFT_ENABLE_CUDA_FORTRAN)
    target_link_options(bench_spline1d_many PRIVATE -cuda -cudaforlibs)
endif()

if(DRAFT_ENABLE_OPENACC)
    set_source_files_properties(src/spline1d_many_openacc.f90 PROPERTIES
        COMPILE_OPTIONS "-acc;-gpu=ccnative;-Minfo=accel"
    )
endif()

if(DRAFT_ENABLE_OMP_TARGET)
    set_source_files_properties(src/spline1d_many_omp_target.f90 PROPERTIES
        COMPILE_OPTIONS "-mp=gpu;-gpu=ccnative;-Minfo=mp"
    )
endif()

if(DRAFT_ENABLE_CUDA_FORTRAN)
    set_source_files_properties(src/spline1d_many_cuda_fortran.f90 PROPERTIES
        COMPILE_OPTIONS "-cuda;-gpu=ccnative"
    )
endif()

if(DRAFT_ENABLE_CUDA_C)
    find_program(NVCC_EXECUTABLE nvcc HINTS /opt/cuda/bin REQUIRED)

    set(CUDA_C_OBJ ${CMAKE_CURRENT_BINARY_DIR}/spline1d_many_cuda.o)
    add_custom_command(
        OUTPUT ${CUDA_C_OBJ}
        COMMAND ${NVCC_EXECUTABLE} -c -O3 -Xcompiler=-fPIC -arch=native
                ${CMAKE_CURRENT_LIST_DIR}/src/spline1d_many_cuda.cu -o ${CUDA_C_OBJ}
        DEPENDS ${CMAKE_CURRENT_LIST_DIR}/src/spline1d_many_cuda.cu
        VERBATIM
    )
    add_library(draft_cuda_c_kernels STATIC ${CUDA_C_OBJ})
    set_target_properties(draft_cuda_c_kernels PROPERTIES LINKER_LANGUAGE C)
    target_link_libraries(draft_spline1d_many PUBLIC draft_cuda_c_kernels)
endif()

if(DRAFT_ENABLE_CUDA_C)
    target_link_libraries(bench_spline1d_many PRIVATE cudart stdc++)
    target_link_directories(bench_spline1d_many PRIVATE /opt/cuda/lib64 /opt/cuda/targets/x86_64-linux/lib)
endif()
