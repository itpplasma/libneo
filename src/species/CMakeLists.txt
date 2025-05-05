add_library(species STATIC
    species.f90
)

set_property(TARGET species PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})

add_library(LIBNEO::species ALIAS species)
target_link_libraries(species neo)
