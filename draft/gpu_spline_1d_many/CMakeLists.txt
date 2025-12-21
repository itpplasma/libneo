cmake_minimum_required(VERSION 3.22)

project(libneo_gpu_spline_1d_many LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2018)

option(DRAFT_ENABLE_OPENACC "Enable OpenACC directives" ON)
set(DRAFT_OPENACC_OFFLOAD "none" CACHE STRING "OpenACC offload target for GNU (none|nvptx)")
set_property(CACHE DRAFT_OPENACC_OFFLOAD PROPERTY STRINGS none nvptx)

set(DRAFT_ROOT ${CMAKE_CURRENT_LIST_DIR})
set(LIBNEO_ROOT ${DRAFT_ROOT}/../..)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(DRAFT_FORTRAN_IS_NVHPC TRUE)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(DRAFT_FORTRAN_IS_GNU TRUE)
else()
    message(FATAL_ERROR "Unsupported Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

set(DRAFT_COMMON_FFLAGS -O3)

set(DRAFT_PREPROCESS_FFLAGS "")
if(DRAFT_FORTRAN_IS_NVHPC)
    list(APPEND DRAFT_PREPROCESS_FFLAGS -Mpreprocess)
elseif(DRAFT_FORTRAN_IS_GNU)
    list(APPEND DRAFT_PREPROCESS_FFLAGS -cpp)
endif()

set(DRAFT_OPENACC_FFLAGS "")
set(DRAFT_OPENACC_LDFLAGS "")
if(DRAFT_ENABLE_OPENACC)
    if(DRAFT_FORTRAN_IS_NVHPC)
        list(APPEND DRAFT_OPENACC_FFLAGS -acc -Minfo=accel)
        list(APPEND DRAFT_OPENACC_LDFLAGS -acc)
        list(APPEND DRAFT_OPENACC_LDFLAGS -gpu=ccnative)
    elseif(DRAFT_FORTRAN_IS_GNU)
        list(APPEND DRAFT_OPENACC_FFLAGS -fopenacc)
        list(APPEND DRAFT_OPENACC_LDFLAGS -fopenacc)
        if(DRAFT_OPENACC_OFFLOAD STREQUAL "nvptx")
            list(APPEND DRAFT_OPENACC_FFLAGS -foffload=nvptx-none)
            list(APPEND DRAFT_OPENACC_LDFLAGS -foffload=nvptx-none)
        elseif(NOT DRAFT_OPENACC_OFFLOAD STREQUAL "none")
            message(FATAL_ERROR "Invalid DRAFT_OPENACC_OFFLOAD='${DRAFT_OPENACC_OFFLOAD}' (expected none|nvptx)")
        endif()
    endif()
endif()

add_library(draft_build_flags INTERFACE)
target_compile_options(draft_build_flags INTERFACE
    $<$<COMPILE_LANGUAGE:Fortran>:${DRAFT_COMMON_FFLAGS}>
    $<$<COMPILE_LANGUAGE:Fortran>:${DRAFT_PREPROCESS_FFLAGS}>
)

if(DRAFT_ENABLE_OPENACC)
    target_compile_definitions(draft_build_flags INTERFACE LIBNEO_ENABLE_OPENACC=1)
    target_compile_options(draft_build_flags INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:${DRAFT_OPENACC_FFLAGS}>)
    target_link_options(draft_build_flags INTERFACE ${DRAFT_OPENACC_LDFLAGS})
endif()

add_library(draft_interpolate STATIC
    ${LIBNEO_ROOT}/src/spl_three_to_five.f90
    ${LIBNEO_ROOT}/src/interpolate/batch_interpolate_types.f90
    ${LIBNEO_ROOT}/src/interpolate/batch_interpolate_1d.f90
    ${LIBNEO_ROOT}/src/interpolate/batch_interpolate_2d.f90
    ${LIBNEO_ROOT}/src/interpolate/batch_interpolate_3d.f90
)
target_include_directories(draft_interpolate PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_link_libraries(draft_interpolate PUBLIC draft_build_flags)

add_library(draft_spline1d_many STATIC
    src/spline1d_many_openacc.f90
    src/util_lcg_rng.f90
)
target_link_libraries(draft_spline1d_many PUBLIC draft_interpolate)

add_executable(bench_spline1d_many src/bench_spline1d_many.f90)
target_link_libraries(bench_spline1d_many PRIVATE draft_spline1d_many)

add_library(draft_spline2d_many STATIC
    src/spline2d_many_openacc.f90
)
target_link_libraries(draft_spline2d_many PUBLIC draft_interpolate)

add_executable(bench_spline2d_many src/bench_spline2d_many.f90)
target_link_libraries(bench_spline2d_many PRIVATE draft_spline2d_many draft_spline1d_many)

add_library(draft_spline3d_many STATIC
    src/spline3d_many_openacc.f90
)
target_link_libraries(draft_spline3d_many PUBLIC draft_interpolate)

add_executable(bench_spline3d_many src/bench_spline3d_many.f90)
target_link_libraries(bench_spline3d_many PRIVATE draft_spline3d_many draft_spline1d_many)

add_library(draft_many_api STATIC
    src/draft_batch_splines_many_api.f90
)
target_link_libraries(draft_many_api PUBLIC draft_spline1d_many draft_spline2d_many draft_spline3d_many)

target_link_libraries(bench_spline1d_many PRIVATE draft_many_api)
target_link_libraries(bench_spline2d_many PRIVATE draft_many_api)
target_link_libraries(bench_spline3d_many PRIVATE draft_many_api)

target_link_libraries(bench_spline1d_many PRIVATE draft_build_flags)
target_link_libraries(bench_spline2d_many PRIVATE draft_build_flags)
target_link_libraries(bench_spline3d_many PRIVATE draft_build_flags)
