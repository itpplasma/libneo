add_library(neo_math STATIC
    bessel_i.f90
    bessel_k.f90
    dawson.f90
    fft.f90
    gauss_kronrod.f90
    gauss_quadrature.f90)

set_property(TARGET neo_math PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})

add_library(LIBNEO::math ALIAS neo_math)

# The math kit relies on IEEE semantics (overflow to Inf, subnormals,
# QUADPACK roundoff detection), which -ffast-math from Release builds breaks.
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(neo_math PRIVATE -fno-fast-math)
endif()

target_link_libraries(neo_math PUBLIC
    neo
    LAPACK::LAPACK
    BLAS::BLAS
)
