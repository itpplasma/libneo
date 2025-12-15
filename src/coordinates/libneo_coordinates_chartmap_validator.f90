submodule(libneo_coordinates) libneo_coordinates_chartmap_validator
    use math_constants, only: TWOPI
    use netcdf, only: NF90_BYTE, NF90_CHAR, NF90_DOUBLE, NF90_GLOBAL, NF90_INT, &
                      NF90_INT64, NF90_MAX_VAR_DIMS, NF90_NOERR, NF90_NOWRITE, &
                      NF90_SHORT, &
                      nf90_close, &
                      nf90_get_att, nf90_get_var, nf90_inq_dimid, nf90_inq_varid, &
                      nf90_inquire_attribute, nf90_inquire_dimension, &
                      nf90_inquire_variable, &
                      nf90_open, nf90_strerror
    implicit none

contains

    module subroutine validate_chartmap_file(filename, ierr, message)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer, parameter :: ok = 0
        integer, parameter :: err_open = 1
        integer, parameter :: err_missing_dim = 2
        integer, parameter :: err_invalid_dim = 3
        integer, parameter :: err_missing_var = 4
        integer, parameter :: err_bad_var = 5
        integer, parameter :: err_missing_units = 6
        integer, parameter :: err_bad_units = 7
        integer, parameter :: err_bad_grid = 8

        integer :: status
        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: len_rho, len_theta, len_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z
        integer :: num_field_periods
        real(dp), allocatable :: rho(:), theta(:), zeta(:)
        real(dp) :: period

        ierr = ok
        message = ""

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            ierr = err_open
            message = "open failed: "//trim(nf90_strerror(status))
            return
        end if

        do
            call require_dim(ncid, "rho", dim_rho, len_rho, ierr, message)
            if (ierr /= ok) exit
            call require_dim(ncid, "theta", dim_theta, len_theta, ierr, message)
            if (ierr /= ok) exit
            call require_dim(ncid, "zeta", dim_zeta, len_zeta, ierr, message)
            if (ierr /= ok) exit

            if (len_rho < 2 .or. len_theta < 2 .or. len_zeta < 1) then
                ierr = err_invalid_dim
                message = "invalid dimension length(s)"
                exit
            end if

            call require_1d_real64(ncid, "rho", dim_rho, var_rho, ierr, message)
            if (ierr /= ok) exit
            call require_1d_real64(ncid, "theta", dim_theta, var_theta, ierr, message)
            if (ierr /= ok) exit
            call require_1d_real64(ncid, "zeta", dim_zeta, var_zeta, ierr, message)
            if (ierr /= ok) exit

            call require_3d_real64(ncid, "x", dim_rho, dim_theta, dim_zeta, var_x, &
                                   ierr, message)
            if (ierr /= ok) exit
            call require_3d_real64(ncid, "y", dim_rho, dim_theta, dim_zeta, var_y, &
                                   ierr, message)
            if (ierr /= ok) exit
            call require_3d_real64(ncid, "z", dim_rho, dim_theta, dim_zeta, var_z, &
                                   ierr, message)
            if (ierr /= ok) exit

            call require_units_cm(ncid, "x", var_x, ierr, message)
            if (ierr /= ok) exit
            call require_units_cm(ncid, "y", var_y, ierr, message)
            if (ierr /= ok) exit
            call require_units_cm(ncid, "z", var_z, ierr, message)
            if (ierr /= ok) exit

            call check_optional_zeta_convention(ncid, ierr, message)
            if (ierr /= ok) exit

            call read_optional_num_field_periods(ncid, num_field_periods, ierr, message)
            if (ierr /= ok) exit

            allocate (rho(len_rho), theta(len_theta), zeta(len_zeta))
            call read_grid(ncid, var_rho, "rho", rho, ierr, message)
            if (ierr /= ok) exit
            call read_grid(ncid, var_theta, "theta", theta, ierr, message)
            if (ierr /= ok) exit
            call read_grid(ncid, var_zeta, "zeta", zeta, ierr, message)
            if (ierr /= ok) exit

            call check_rho_grid(rho, ierr, message)
            if (ierr /= ok) exit

            call check_periodic_grid("theta", theta, TWOPI, ierr, message)
            if (ierr /= ok) exit

            period = TWOPI/real(num_field_periods, dp)
            call check_periodic_grid("zeta", zeta, period, ierr, message)
            if (ierr /= ok) exit

            exit
        end do

        status = nf90_close(ncid)
    end subroutine validate_chartmap_file

    subroutine check_optional_zeta_convention(ncid, ierr, message)
        integer, intent(in) :: ncid
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status
        integer :: att_type
        integer :: att_len
        character(len=:), allocatable :: value

        ierr = 0
        message = ""

        status = nf90_inquire_attribute(ncid, NF90_GLOBAL, "zeta_convention", &
                                        xtype=att_type, len=att_len)
        if (status /= NF90_NOERR) return
        if (att_type /= NF90_CHAR .or. att_len < 1) then
            ierr = 5
            message = "zeta_convention must be a global string attribute"
            return
        end if

        allocate (character(len=att_len) :: value)
        status = nf90_get_att(ncid, NF90_GLOBAL, "zeta_convention", value)
        if (status /= NF90_NOERR) then
            ierr = 5
            message = "could not read zeta_convention"
            return
        end if

        select case (trim(value))
        case ("cyl", "vmec", "boozer", "unknown")
        case default
            ierr = 5
            message = "invalid zeta_convention"
        end select
    end subroutine check_optional_zeta_convention

    subroutine require_dim(ncid, name, dimid, dimlen, ierr, message)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(out) :: dimid
        integer, intent(out) :: dimlen
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status

        ierr = 0
        message = ""

        status = nf90_inq_dimid(ncid, trim(name), dimid)
        if (status /= NF90_NOERR) then
            ierr = 2
            message = "missing dimension "//trim(name)
            return
        end if

        status = nf90_inquire_dimension(ncid, dimid, len=dimlen)
        if (status /= NF90_NOERR) then
            ierr = 3
            message = "invalid dimension "//trim(name)
            return
        end if
    end subroutine require_dim

    subroutine require_1d_real64(ncid, name, dimid, varid, ierr, message)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(in) :: dimid
        integer, intent(out) :: varid
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status, xtype, ndims
        integer :: dimids(NF90_MAX_VAR_DIMS)

        ierr = 0
        message = ""

        status = nf90_inq_varid(ncid, trim(name), varid)
        if (status /= NF90_NOERR) then
            ierr = 4
            message = "missing variable "//trim(name)
            return
        end if

        status = nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=ndims, &
                                       dimids=dimids)
        if (status /= NF90_NOERR) then
            ierr = 5
            message = "could not inquire variable "//trim(name)
            return
        end if

        if (xtype /= NF90_DOUBLE .or. ndims /= 1 .or. dimids(1) /= dimid) then
            ierr = 5
            message = "variable "//trim(name)//" must be float64 with correct dimension"
            return
        end if
    end subroutine require_1d_real64

    subroutine require_3d_real64(ncid, name, dim_rho, dim_theta, dim_zeta, varid, &
                                 ierr, message)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(in) :: dim_rho, dim_theta, dim_zeta
        integer, intent(out) :: varid
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status, xtype, ndims
        integer :: dimids(NF90_MAX_VAR_DIMS)

        ierr = 0
        message = ""

        status = nf90_inq_varid(ncid, trim(name), varid)
        if (status /= NF90_NOERR) then
            ierr = 4
            message = "missing variable "//trim(name)
            return
        end if

        status = nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=ndims, &
                                       dimids=dimids)
        if (status /= NF90_NOERR) then
            ierr = 5
            message = "could not inquire variable "//trim(name)
            return
        end if

        if (xtype /= NF90_DOUBLE .or. ndims /= 3) then
            ierr = 5
            message = "variable "//trim(name)//" must be float64 with 3 dimensions"
            return
        end if

        if (dimids(1) /= dim_rho .or. dimids(2) /= dim_theta .or. dimids(3) /= &
            dim_zeta) then
            ierr = 5
            message = "variable "//trim(name)// &
                      " has wrong dimension order; expected (rho,theta,zeta)"
            return
        end if
    end subroutine require_3d_real64

    subroutine require_units_cm(ncid, name, varid, ierr, message)
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        integer, intent(in) :: varid
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status, att_len, att_type
        character(len=:), allocatable :: units

        ierr = 0
        message = ""

        status = nf90_inquire_attribute(ncid, varid, "units", xtype=att_type, &
                                        len=att_len)
        if (status /= NF90_NOERR) then
            ierr = 6
            message = "missing units attribute for "//trim(name)
            return
        end if

        if (att_type /= NF90_CHAR .or. att_len < 1) then
            ierr = 7
            message = "invalid units attribute for "//trim(name)
            return
        end if

        allocate (character(len=att_len) :: units)
        status = nf90_get_att(ncid, varid, "units", units)
        if (status /= NF90_NOERR) then
            ierr = 7
            message = "could not read units attribute for "//trim(name)
            return
        end if

        if (trim(units) /= "cm") then
            ierr = 7
            message = "units for "//trim(name)//" must be cm"
            return
        end if
    end subroutine require_units_cm

    subroutine read_optional_num_field_periods(ncid, num_field_periods, ierr, message)
        integer, intent(in) :: ncid
        integer, intent(out) :: num_field_periods
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status, varid, xtype, ndims
        character(len=*), parameter :: var_primary = "num_field_periods"
        character(len=*), parameter :: var_legacy = "nfp"

        ierr = 0
        message = ""
        num_field_periods = 1

        status = nf90_inq_varid(ncid, var_primary, varid)
        if (status /= NF90_NOERR) then
            status = nf90_inq_varid(ncid, var_legacy, varid)
            if (status /= NF90_NOERR) return
        end if

        status = nf90_inquire_variable(ncid, varid, xtype=xtype, ndims=ndims)
        if (status /= NF90_NOERR) then
            ierr = 5
            message = "could not inquire variable num_field_periods"
            return
        end if

        if (ndims /= 0 .or. .not. is_integer_xtype(xtype)) then
            ierr = 5
            message = "num_field_periods must be a scalar integer variable"
            return
        end if

        status = nf90_get_var(ncid, varid, num_field_periods)
        if (status /= NF90_NOERR) then
            ierr = 5
            message = "could not read num_field_periods"
            return
        end if

        if (num_field_periods < 1) then
            ierr = 5
            message = "num_field_periods must be >= 1"
            return
        end if
    end subroutine read_optional_num_field_periods

    subroutine read_grid(ncid, varid, name, x, ierr, message)
        integer, intent(in) :: ncid
        integer, intent(in) :: varid
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: status

        ierr = 0
        message = ""

        status = nf90_get_var(ncid, varid, x)
        if (status /= NF90_NOERR) then
            ierr = 5
            message = "could not read "//trim(name)
        end if
    end subroutine read_grid

    subroutine check_rho_grid(rho, ierr, message)
        real(dp), intent(in) :: rho(:)
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        real(dp), parameter :: tol = 1.0e-12_dp

        ierr = 0
        message = ""

        if (.not. is_strictly_increasing(rho)) then
            ierr = 8
            message = "rho grid must be strictly increasing"
            return
        end if

        if (abs(rho(1) - 0.0_dp) > tol .or. abs(rho(size(rho)) - 1.0_dp) > tol) then
            ierr = 8
            message = "rho must span [0,1] including endpoints"
            return
        end if
    end subroutine check_rho_grid

    subroutine check_periodic_grid(name, x, period, ierr, message)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: period
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        real(dp), parameter :: tol = 1.0e-12_dp
        real(dp) :: step_expected

        ierr = 0
        message = ""

        if (size(x) == 1) then
            if (abs(x(1) - 0.0_dp) > tol) then
                ierr = 8
                message = trim(name)//" must start at 0"
            end if
            return
        end if

        if (.not. is_strictly_increasing(x)) then
            ierr = 8
            message = trim(name)//" grid must be strictly increasing"
            return
        end if

        step_expected = period/real(size(x), dp)
        if (abs(x(1) - 0.0_dp) > tol) then
            ierr = 8
            message = trim(name)//" must start at 0"
            return
        end if

        if (.not. is_uniform_step(x, step_expected, tol)) then
            ierr = 8
            message = trim(name)//" must be uniform with endpoint excluded"
            return
        end if

        if (x(size(x)) > (period - step_expected + tol)) then
            ierr = 8
            message = trim(name)//" must exclude the period endpoint"
            return
        end if
    end subroutine check_periodic_grid

    logical function is_integer_xtype(xtype)
        integer, intent(in) :: xtype
        is_integer_xtype = (xtype == NF90_BYTE) .or. (xtype == NF90_SHORT) .or. &
                           (xtype == NF90_INT) .or. (xtype == NF90_INT64)
    end function is_integer_xtype

    logical function is_strictly_increasing(x)
        real(dp), intent(in) :: x(:)
        integer :: i

        is_strictly_increasing = .true.
        do i = 2, size(x)
            if (x(i) <= x(i - 1)) then
                is_strictly_increasing = .false.
                return
            end if
        end do
    end function is_strictly_increasing

    logical function is_uniform_step(x, step, tol)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: step
        real(dp), intent(in) :: tol

        real(dp) :: dx_max

        if (size(x) < 2) then
            is_uniform_step = .true.
            return
        end if

        dx_max = maxval(abs((x(2:) - x(:size(x) - 1)) - step))
        is_uniform_step = (dx_max <= 10.0_dp*tol)
    end function is_uniform_step

end submodule libneo_coordinates_chartmap_validator
