cmake_minimum_required(VERSION 3.22)

project(libneo_bench_spline_many LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2018)

option(BENCH_ENABLE_OPENACC "Enable OpenACC directives" ON)

set(BENCH_ROOT ${CMAKE_CURRENT_LIST_DIR})
set(LIBNEO_ROOT ${BENCH_ROOT}/../..)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(BENCH_FORTRAN_IS_NVHPC TRUE)
else()
    message(FATAL_ERROR "Benchmarks support NVHPC nvfortran only (got: ${CMAKE_Fortran_COMPILER_ID})")
endif()

set(BENCH_COMMON_FFLAGS -O3)

set(BENCH_PREPROCESS_FFLAGS "")
if(BENCH_FORTRAN_IS_NVHPC)
    list(APPEND BENCH_PREPROCESS_FFLAGS -Mpreprocess)
endif()

set(BENCH_OPENACC_FFLAGS "")
set(BENCH_OPENACC_LDFLAGS "")
if(BENCH_ENABLE_OPENACC)
    if(BENCH_FORTRAN_IS_NVHPC)
        list(APPEND BENCH_OPENACC_FFLAGS -acc -Minfo=accel)
        list(APPEND BENCH_OPENACC_LDFLAGS -acc)
        list(APPEND BENCH_OPENACC_LDFLAGS -gpu=ccnative)
    endif()
endif()

add_library(bench_build_flags INTERFACE)
target_compile_options(bench_build_flags INTERFACE
    $<$<COMPILE_LANGUAGE:Fortran>:${BENCH_COMMON_FFLAGS}>
    $<$<COMPILE_LANGUAGE:Fortran>:${BENCH_PREPROCESS_FFLAGS}>
)

if(BENCH_ENABLE_OPENACC)
    target_compile_definitions(bench_build_flags INTERFACE LIBNEO_ENABLE_OPENACC=1)
    target_compile_options(bench_build_flags INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:${BENCH_OPENACC_FFLAGS}>)
    target_link_options(bench_build_flags INTERFACE ${BENCH_OPENACC_LDFLAGS})
endif()

add_library(bench_interpolate STATIC
    ${LIBNEO_ROOT}/src/spl_three_to_five.f90
    ${LIBNEO_ROOT}/src/interpolate/acc_spline_build_order5.f90
    ${LIBNEO_ROOT}/src/interpolate/batch_interpolate_types.f90
    ${LIBNEO_ROOT}/src/interpolate/batch_interpolate_1d.f90
    ${LIBNEO_ROOT}/src/interpolate/batch_interpolate_2d.f90
    ${LIBNEO_ROOT}/src/interpolate/batch_interpolate_3d.f90
)
target_include_directories(bench_interpolate PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_link_libraries(bench_interpolate PUBLIC bench_build_flags)

add_library(bench_spline1d_many_lib STATIC
    src/spline1d_many_offload.f90
    src/util_lcg_rng.f90
    src/util_timer.f90
    src/util_bench_args.f90
)
target_link_libraries(bench_spline1d_many_lib PUBLIC bench_interpolate)

add_executable(bench_spline1d_many src/bench_spline1d_many.f90)
target_link_libraries(bench_spline1d_many PRIVATE bench_spline1d_many_lib)

add_library(bench_spline2d_many_lib STATIC
    src/spline2d_many_offload.f90
)
target_link_libraries(bench_spline2d_many_lib PUBLIC bench_interpolate)

add_executable(bench_spline2d_many src/bench_spline2d_many.f90)
target_link_libraries(bench_spline2d_many PRIVATE bench_spline2d_many_lib bench_spline1d_many_lib)

add_library(bench_spline3d_many_lib STATIC
    src/spline3d_many_offload.f90
)
target_link_libraries(bench_spline3d_many_lib PUBLIC bench_interpolate)

add_executable(bench_spline3d_many src/bench_spline3d_many.f90)
target_link_libraries(bench_spline3d_many PRIVATE bench_spline3d_many_lib bench_spline1d_many_lib)

add_library(bench_many_api STATIC
    src/batch_splines_many_api.f90
)
target_link_libraries(bench_many_api PUBLIC bench_spline1d_many_lib bench_spline2d_many_lib bench_spline3d_many_lib)

target_link_libraries(bench_spline1d_many PRIVATE bench_many_api)
target_link_libraries(bench_spline2d_many PRIVATE bench_many_api)
target_link_libraries(bench_spline3d_many PRIVATE bench_many_api)

target_link_libraries(bench_spline1d_many PRIVATE bench_build_flags)
target_link_libraries(bench_spline2d_many PRIVATE bench_build_flags)
target_link_libraries(bench_spline3d_many PRIVATE bench_build_flags)
