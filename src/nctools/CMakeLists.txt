add_library(nctools STATIC
    nctools_module.f90
)

# Link NetCDF - prefer imported targets, fall back to nf-config
if(TARGET netcdf::netcdff)
    target_link_libraries(nctools PUBLIC netcdf::netcdff)
elseif(TARGET netcdf::netcdf)
    target_link_libraries(nctools PUBLIC netcdf::netcdf netcdff)
elseif(NETCDF_FLIBS_LIST)
    # System NetCDF found by DetectDependencies.cmake
    target_include_directories(nctools PUBLIC ${NETCDFINCLUDE_DIR})
    if(NETCDFLIB_DIR)
        target_link_directories(nctools PUBLIC ${NETCDFLIB_DIR})
    endif()
    # Also add nc-config libdir for the C netcdf library
    find_program(_nctools_nc_config nc-config)
    if(_nctools_nc_config)
        execute_process(COMMAND nc-config --libdir
            OUTPUT_VARIABLE _nctools_nc_libdir OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(_nctools_nc_libdir)
            target_link_directories(nctools PUBLIC ${_nctools_nc_libdir})
        endif()
    endif()
    target_link_libraries(nctools PUBLIC ${NETCDF_FLIBS_LIST})
else()
    # Last resort: probe nf-config/nc-config ourselves
    find_program(_nctools_nf_config nf-config)
    if(_nctools_nf_config)
        execute_process(COMMAND nf-config --includedir
            OUTPUT_VARIABLE _nctools_nc_inc OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND nf-config --flibs
            OUTPUT_VARIABLE _nctools_nc_flibs OUTPUT_STRIP_TRAILING_WHITESPACE)
        string(REPLACE " " ";" _nctools_nc_flibs_list "${_nctools_nc_flibs}")
        target_include_directories(nctools PUBLIC ${_nctools_nc_inc})
        # Also add nc-config libdir for the C library (-lnetcdf)
        find_program(_nctools_nc_config nc-config)
        if(_nctools_nc_config)
            execute_process(COMMAND nc-config --libdir
                OUTPUT_VARIABLE _nctools_nc_libdir OUTPUT_STRIP_TRAILING_WHITESPACE)
            target_link_directories(nctools PUBLIC ${_nctools_nc_libdir})
        endif()
        target_link_libraries(nctools PUBLIC ${_nctools_nc_flibs_list})
    else()
        target_link_libraries(nctools PUBLIC netcdf netcdff)
    endif()
endif()

set_property(TARGET nctools PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::nctools ALIAS nctools)
