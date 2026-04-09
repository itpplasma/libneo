add_library(nctools STATIC
    nctools_module.f90
)

# Link NetCDF
if(TARGET netcdf::netcdff)
    target_link_libraries(nctools PUBLIC netcdf::netcdff)
elseif(TARGET netcdf::netcdf)
    target_link_libraries(nctools PUBLIC netcdf::netcdf netcdff)
else()
    target_link_libraries(nctools PUBLIC netcdf netcdff)
endif()

set_property(TARGET nctools PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::nctools ALIAS nctools)
