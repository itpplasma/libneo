add_library(coordinates STATIC
    cylindrical_cartesian.f90
    vmec_coordinates.f90
    geoflux_coordinates.f90
    ../magfie/geqdsk_tools.f90
    libneo_coordinate_conventions.f90
    libneo_coordinates_base.f90
    libneo_coordinates_vmec.f90
    libneo_coordinates_geoflux.f90
    libneo_coordinates_chartmap.f90
    libneo_coordinates_validator.f90
    libneo_coordinates_file_detection.f90
    libneo_coordinates.f90
    libneo_chartmap_vmec_generator.f90
    libneo_coordinates_mapping.f90
)

target_link_libraries(coordinates PUBLIC
    neo_util
    nctools
    vmec_support
    interpolate
)

# Link NetCDF directly (chartmap uses netcdf module)
if(TARGET netcdf::netcdff)
    target_link_libraries(coordinates PRIVATE netcdf::netcdff)
elseif(TARGET netcdf::netcdf)
    target_link_libraries(coordinates PRIVATE netcdf::netcdf netcdff)
else()
    target_link_libraries(coordinates PRIVATE netcdf netcdff)
endif()

set_property(TARGET coordinates PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::coordinates ALIAS coordinates)
