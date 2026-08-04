add_library(nctools STATIC
    nctools_module.f90
)

target_link_libraries(nctools PUBLIC fortio::fortio)

set_property(TARGET nctools PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::nctools ALIAS nctools)
