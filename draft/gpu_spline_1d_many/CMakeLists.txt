cmake_minimum_required(VERSION 3.22)

project(libneo_gpu_spline_1d_many LANGUAGES C Fortran)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_C_STANDARD 11)

option(DRAFT_ENABLE_OPENACC "Build OpenACC variant (NVHPC)" ON)
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
    ${LIBNEO_ROOT}/src/interpolate/batch_interpolate_2d.f90
    ${LIBNEO_ROOT}/src/interpolate/batch_interpolate_3d.f90
)
target_include_directories(draft_interpolate PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_compile_options(draft_interpolate PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-O3>)

add_library(draft_spline1d_many STATIC
    src/spline1d_many_cpu.f90
    src/spline1d_many_openacc.f90
    src/spline1d_many_cuda_fortran.cuf
    src/spline1d_many_cuda_c_wrapper.f90
    src/util_lcg_rng.f90
)
target_link_libraries(draft_spline1d_many PUBLIC draft_interpolate)
target_compile_options(draft_spline1d_many PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-O3>)

add_executable(bench_spline1d_many src/bench_spline1d_many.f90)
target_link_libraries(bench_spline1d_many PRIVATE draft_spline1d_many)
target_compile_options(bench_spline1d_many PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-O3>)

target_link_options(bench_spline1d_many PRIVATE
    -gpu=ccnative
)

if(DRAFT_ENABLE_OPENACC)
    target_link_options(bench_spline1d_many PRIVATE -acc)
endif()

if(DRAFT_ENABLE_CUDA_FORTRAN)
    target_link_options(bench_spline1d_many PRIVATE -cuda -cudaforlibs)
endif()

if(DRAFT_ENABLE_OPENACC)
    set_source_files_properties(src/spline1d_many_openacc.f90 PROPERTIES
        COMPILE_DEFINITIONS "SPLINE1D_BACKEND_OPENACC"
        COMPILE_OPTIONS "-acc;-gpu=ccnative;-Minfo=accel"
    )
endif()

if(DRAFT_ENABLE_CUDA_FORTRAN)
    set_source_files_properties(src/spline1d_many_cuda_fortran.cuf PROPERTIES
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

add_library(draft_spline2d_many STATIC
    src/spline2d_many_cpu.f90
    src/spline2d_many_openacc.f90
    src/spline2d_many_cuda_fortran.cuf
    src/spline2d_many_cuda_c_wrapper.f90
)
target_link_libraries(draft_spline2d_many PUBLIC draft_interpolate)
target_compile_options(draft_spline2d_many PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-O3>)

add_executable(bench_spline2d_many src/bench_spline2d_many.f90)
target_link_libraries(bench_spline2d_many PRIVATE draft_spline2d_many draft_spline1d_many)
target_compile_options(bench_spline2d_many PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-O3>)
target_link_options(bench_spline2d_many PRIVATE -gpu=ccnative)

if(DRAFT_ENABLE_OPENACC)
    target_link_options(bench_spline2d_many PRIVATE -acc)
endif()

if(DRAFT_ENABLE_CUDA_FORTRAN)
    target_link_options(bench_spline2d_many PRIVATE -cuda -cudaforlibs)
endif()

if(DRAFT_ENABLE_OPENACC)
    set_source_files_properties(src/spline2d_many_openacc.f90 PROPERTIES
        COMPILE_DEFINITIONS "SPLINE2D_BACKEND_OPENACC"
        COMPILE_OPTIONS "-acc;-gpu=ccnative;-Minfo=accel"
    )
endif()

if(DRAFT_ENABLE_CUDA_FORTRAN)
    set_source_files_properties(src/spline2d_many_cuda_fortran.cuf PROPERTIES
        COMPILE_OPTIONS "-cuda;-gpu=ccnative"
    )
endif()

if(DRAFT_ENABLE_CUDA_C)
    set(CUDA2D_C_OBJ ${CMAKE_CURRENT_BINARY_DIR}/spline2d_many_cuda.o)
    add_custom_command(
        OUTPUT ${CUDA2D_C_OBJ}
        COMMAND ${NVCC_EXECUTABLE} -c -O3 -Xcompiler=-fPIC -arch=native
                ${CMAKE_CURRENT_LIST_DIR}/src/spline2d_many_cuda.cu -o ${CUDA2D_C_OBJ}
        DEPENDS ${CMAKE_CURRENT_LIST_DIR}/src/spline2d_many_cuda.cu
        VERBATIM
    )
    add_library(draft_cuda2d_c_kernels STATIC ${CUDA2D_C_OBJ})
    set_target_properties(draft_cuda2d_c_kernels PROPERTIES LINKER_LANGUAGE C)
    target_link_libraries(draft_spline2d_many PUBLIC draft_cuda2d_c_kernels)
    target_link_libraries(bench_spline2d_many PRIVATE cudart stdc++)
    target_link_directories(bench_spline2d_many PRIVATE /opt/cuda/lib64 /opt/cuda/targets/x86_64-linux/lib)
endif()

add_library(draft_spline3d_many STATIC
    src/spline3d_many_cpu.f90
    src/spline3d_many_openacc.f90
    src/spline3d_many_cuda_fortran.cuf
    src/spline3d_many_cuda_c_wrapper.f90
)
target_link_libraries(draft_spline3d_many PUBLIC draft_interpolate)
target_compile_options(draft_spline3d_many PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-O3>)

add_executable(bench_spline3d_many src/bench_spline3d_many.f90)
target_link_libraries(bench_spline3d_many PRIVATE draft_spline3d_many draft_spline1d_many)
target_compile_options(bench_spline3d_many PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-O3>)
target_link_options(bench_spline3d_many PRIVATE -gpu=ccnative)

if(DRAFT_ENABLE_OPENACC)
    target_link_options(bench_spline3d_many PRIVATE -acc)
    set_source_files_properties(src/spline3d_many_openacc.f90 PROPERTIES
        COMPILE_DEFINITIONS "SPLINE3D_BACKEND_OPENACC"
        COMPILE_OPTIONS "-acc;-gpu=ccnative;-Minfo=accel"
    )
endif()

if(DRAFT_ENABLE_CUDA_FORTRAN)
    target_link_options(bench_spline3d_many PRIVATE -cuda -cudaforlibs)
    set_source_files_properties(src/spline3d_many_cuda_fortran.cuf PROPERTIES
        COMPILE_OPTIONS "-cuda;-gpu=ccnative"
    )
endif()

if(DRAFT_ENABLE_CUDA_C)
    set(CUDA3D_C_OBJ ${CMAKE_CURRENT_BINARY_DIR}/spline3d_many_cuda.o)
    add_custom_command(
        OUTPUT ${CUDA3D_C_OBJ}
        COMMAND ${NVCC_EXECUTABLE} -c -O3 -Xcompiler=-fPIC -arch=native
                ${CMAKE_CURRENT_LIST_DIR}/src/spline3d_many_cuda.cu -o ${CUDA3D_C_OBJ}
        DEPENDS ${CMAKE_CURRENT_LIST_DIR}/src/spline3d_many_cuda.cu
        VERBATIM
    )
    add_library(draft_cuda3d_c_kernels STATIC ${CUDA3D_C_OBJ})
    set_target_properties(draft_cuda3d_c_kernels PROPERTIES LINKER_LANGUAGE C)
    target_link_libraries(draft_spline3d_many PUBLIC draft_cuda3d_c_kernels)
    target_link_libraries(bench_spline3d_many PRIVATE cudart stdc++)
    target_link_directories(bench_spline3d_many PRIVATE /opt/cuda/lib64 /opt/cuda/targets/x86_64-linux/lib)
endif()

add_library(draft_many_api STATIC
    src/draft_batch_splines_many_api.f90
)
target_link_libraries(draft_many_api PUBLIC draft_spline1d_many draft_spline2d_many draft_spline3d_many)
target_compile_options(draft_many_api PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-O3>)

target_link_libraries(bench_spline1d_many PRIVATE draft_many_api)
target_link_libraries(bench_spline2d_many PRIVATE draft_many_api)
target_link_libraries(bench_spline3d_many PRIVATE draft_many_api)
