cmake_minimum_required(VERSION 3.22)

project(libneo_gpu_spline_1d_many LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2018)

option(DRAFT_ENABLE_OPENACC "Build OpenACC variant (NVHPC)" ON)

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
    src/spline1d_many_openacc.f90
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

if(DRAFT_ENABLE_OPENACC)
    set_source_files_properties(src/spline1d_many_openacc.f90 PROPERTIES
        COMPILE_OPTIONS "-acc;-gpu=ccnative;-Minfo=accel"
    )
endif()

add_library(draft_spline2d_many STATIC
    src/spline2d_many_openacc.f90
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

if(DRAFT_ENABLE_OPENACC)
    set_source_files_properties(src/spline2d_many_openacc.f90 PROPERTIES
        COMPILE_OPTIONS "-acc;-gpu=ccnative;-Minfo=accel"
    )
endif()

add_library(draft_spline3d_many STATIC
    src/spline3d_many_openacc.f90
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
        COMPILE_OPTIONS "-acc;-gpu=ccnative;-Minfo=accel"
    )
endif()

add_library(draft_many_api STATIC
    src/draft_batch_splines_many_api.f90
)
target_link_libraries(draft_many_api PUBLIC draft_spline1d_many draft_spline2d_many draft_spline3d_many)
target_compile_options(draft_many_api PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-O3>)

target_link_libraries(bench_spline1d_many PRIVATE draft_many_api)
target_link_libraries(bench_spline2d_many PRIVATE draft_many_api)
target_link_libraries(bench_spline3d_many PRIVATE draft_many_api)
