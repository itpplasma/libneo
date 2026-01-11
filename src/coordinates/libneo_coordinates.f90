module libneo_coordinates
    use libneo_coordinates_base
    use libneo_coordinates_vmec
    use libneo_coordinates_geoflux
    use libneo_coordinates_chartmap
    use libneo_coordinates_validator
    use libneo_coordinates_file_detection
    implicit none

    public :: coordinate_system_t
    public :: UNKNOWN, CYL, VMEC, BOOZER
    public :: RHO_TOR, RHO_POL, PSI_TOR_NORM, PSI_POL_NORM

    public :: chartmap_netcdf_spec, chartmap_from_cyl_ierr_spec
    public :: chartmap_from_cyl_ok, chartmap_from_cyl_err_max_iter
    public :: chartmap_from_cyl_err_singular, chartmap_from_cyl_err_out_of_bounds
    public :: chartmap_from_cyl_err_invalid
    public :: refcoords_file_unknown, refcoords_file_chartmap, refcoords_file_vmec_wout

    public :: vmec_coordinate_system_t, make_vmec_coordinate_system
    public :: geoflux_coordinate_system_t, make_geoflux_coordinate_system
    public :: chartmap_coordinate_system_t, make_chartmap_coordinate_system

    public :: validate_chartmap_file
    public :: detect_refcoords_file_type

end module libneo_coordinates
