program test_chartmap_validator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: validate_chartmap_file
    use math_constants, only: TWOPI
    use netcdf
    implicit none

    integer :: ierr
    integer :: nerrors
    character(len=2048) :: message

    nerrors = 0

    call write_good_file("chartmap_good.nc")
    call validate_chartmap_file("chartmap_good.nc", ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: validator rejected a known-good chartmap file"
        print *, trim(message)
        nerrors = nerrors + 1
    else
        print *, "  PASS: validator accepted a known-good chartmap file"
    end if

    call write_missing_x_file("chartmap_missing_x.nc")
    call validate_chartmap_file("chartmap_missing_x.nc", ierr, message)
    if (ierr == 0 .or. index(message, "missing variable x") == 0) then
        print *, "  FAIL: expected missing-x error, got code ", ierr
        print *, trim(message)
        nerrors = nerrors + 1
    else
        print *, "  PASS: missing-x file rejected with expected message"
    end if

    call write_missing_units_file("chartmap_missing_units.nc")
    call validate_chartmap_file("chartmap_missing_units.nc", ierr, message)
    if (ierr == 0 .or. index(message, "units") == 0) then
        print *, "  FAIL: expected missing-units error, got code ", ierr
        print *, trim(message)
        nerrors = nerrors + 1
    else
        print *, "  PASS: missing-units file rejected with expected message"
    end if

    call write_wrong_dims_file("chartmap_wrong_dims.nc")
    call validate_chartmap_file("chartmap_wrong_dims.nc", ierr, message)
    if (ierr == 0 .or. index(message, "dimension order") == 0) then
        print *, "  FAIL: expected wrong-dims error, got code ", ierr
        print *, trim(message)
        nerrors = nerrors + 1
    else
        print *, "  PASS: wrong-dims file rejected with expected message"
    end if

    call write_zeta_endpoint_file("chartmap_zeta_endpoint.nc")
    call validate_chartmap_file("chartmap_zeta_endpoint.nc", ierr, message)
    if (ierr == 0 .or. index(message, "zeta") == 0) then
        print *, "  FAIL: expected zeta-endpoint error, got code ", ierr
        print *, trim(message)
        nerrors = nerrors + 1
    else
        print *, "  PASS: zeta-endpoint file rejected with expected message"
    end if

    call write_unknown_convention_file("chartmap_unknown_convention.nc")
    call validate_chartmap_file("chartmap_unknown_convention.nc", ierr, message)
    if (ierr == 0 .or. index(message, "zeta_convention") == 0) then
        print *, "  FAIL: expected unknown zeta_convention error, got code ", ierr
        print *, trim(message)
        nerrors = nerrors + 1
    else
        print *, "  PASS: unknown zeta_convention rejected with expected message"
    end if

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) detected in chartmap validator tests"
        error stop 1
    end if

contains

    subroutine nc_check(status)
        integer, intent(in) :: status
        if (status /= NF90_NOERR) then
            print *, trim(nf90_strerror(status))
            error stop 2
        end if
    end subroutine nc_check

    subroutine define_grids(nrho, ntheta, nzeta, rho, theta, zeta)
        integer, intent(in) :: nrho, ntheta, nzeta
        real(dp), intent(out) :: rho(nrho), theta(ntheta), zeta(nzeta)
        integer :: i

        do i = 1, nrho
            rho(i) = real(i - 1, dp)/real(nrho - 1, dp)
        end do
        do i = 1, ntheta
            theta(i) = TWOPI*real(i - 1, dp)/real(ntheta, dp)
        end do
        do i = 1, nzeta
            zeta(i) = TWOPI*real(i - 1, dp)/real(nzeta, dp)
        end do
    end subroutine define_grids

    subroutine write_missing_x_file(filename)
        character(len=*), intent(in) :: filename
        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_y, var_z, var_num_field_periods
        integer, parameter :: nrho = 3, ntheta = 4, nzeta = 5
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp) :: y(nrho, ntheta, nzeta), z(nrho, ntheta, nzeta)

        call define_grids(nrho, ntheta, nzeta, rho, theta, zeta)
        y = 0.0_dp
        z = 0.0_dp

        call nc_check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "zeta_convention", "cyl"))
        call nc_check(nf90_def_dim(ncid, "rho", nrho, dim_rho))
        call nc_check(nf90_def_dim(ncid, "theta", ntheta, dim_theta))
        call nc_check(nf90_def_dim(ncid, "zeta", nzeta, dim_zeta))
        call nc_check(nf90_def_var(ncid, "rho", NF90_DOUBLE, [dim_rho], var_rho))
        call nc_check(nf90_def_var(ncid, "theta", NF90_DOUBLE, [dim_theta], var_theta))
        call nc_check(nf90_def_var(ncid, "zeta", NF90_DOUBLE, [dim_zeta], var_zeta))
        call nc_check(nf90_def_var(ncid, "y", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_y))
        call nc_check(nf90_def_var(ncid, "z", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_z))
        call nc_check(nf90_def_var(ncid, "num_field_periods", NF90_INT, &
                                   var_num_field_periods))
        call nc_check(nf90_put_att(ncid, var_y, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_z, "units", "cm"))
        call nc_check(nf90_enddef(ncid))

        call nc_check(nf90_put_var(ncid, var_rho, rho))
        call nc_check(nf90_put_var(ncid, var_theta, theta))
        call nc_check(nf90_put_var(ncid, var_zeta, zeta))
        call nc_check(nf90_put_var(ncid, var_y, y))
        call nc_check(nf90_put_var(ncid, var_z, z))
        call nc_check(nf90_put_var(ncid, var_num_field_periods, 1))
        call nc_check(nf90_close(ncid))
    end subroutine write_missing_x_file

    subroutine write_good_file(filename)
        character(len=*), intent(in) :: filename
        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z, &
                   var_num_field_periods
        integer, parameter :: nrho = 3, ntheta = 4, nzeta = 5
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), z(nrho, &
                                                                      ntheta, nzeta)

        call define_grids(nrho, ntheta, nzeta, rho, theta, zeta)
        x = 0.0_dp
        y = 0.0_dp
        z = 0.0_dp

        call nc_check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "zeta_convention", "cyl"))
        call nc_check(nf90_def_dim(ncid, "rho", nrho, dim_rho))
        call nc_check(nf90_def_dim(ncid, "theta", ntheta, dim_theta))
        call nc_check(nf90_def_dim(ncid, "zeta", nzeta, dim_zeta))
        call nc_check(nf90_def_var(ncid, "rho", NF90_DOUBLE, [dim_rho], var_rho))
        call nc_check(nf90_def_var(ncid, "theta", NF90_DOUBLE, [dim_theta], var_theta))
        call nc_check(nf90_def_var(ncid, "zeta", NF90_DOUBLE, [dim_zeta], var_zeta))
        call nc_check(nf90_def_var(ncid, "x", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_x))
        call nc_check(nf90_def_var(ncid, "y", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_y))
        call nc_check(nf90_def_var(ncid, "z", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_z))
        call nc_check(nf90_def_var(ncid, "num_field_periods", NF90_INT, &
                                   var_num_field_periods))
        call nc_check(nf90_put_att(ncid, var_x, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_y, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_z, "units", "cm"))
        call nc_check(nf90_enddef(ncid))

        call nc_check(nf90_put_var(ncid, var_rho, rho))
        call nc_check(nf90_put_var(ncid, var_theta, theta))
        call nc_check(nf90_put_var(ncid, var_zeta, zeta))
        call nc_check(nf90_put_var(ncid, var_x, x))
        call nc_check(nf90_put_var(ncid, var_y, y))
        call nc_check(nf90_put_var(ncid, var_z, z))
        call nc_check(nf90_put_var(ncid, var_num_field_periods, 1))
        call nc_check(nf90_close(ncid))
    end subroutine write_good_file

    subroutine write_missing_units_file(filename)
        character(len=*), intent(in) :: filename
        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z
        integer, parameter :: nrho = 3, ntheta = 4, nzeta = 5
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), z(nrho, &
                                                                      ntheta, nzeta)

        call define_grids(nrho, ntheta, nzeta, rho, theta, zeta)
        x = 0.0_dp
        y = 0.0_dp
        z = 0.0_dp

        call nc_check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "zeta_convention", "cyl"))
        call nc_check(nf90_def_dim(ncid, "rho", nrho, dim_rho))
        call nc_check(nf90_def_dim(ncid, "theta", ntheta, dim_theta))
        call nc_check(nf90_def_dim(ncid, "zeta", nzeta, dim_zeta))
        call nc_check(nf90_def_var(ncid, "rho", NF90_DOUBLE, [dim_rho], var_rho))
        call nc_check(nf90_def_var(ncid, "theta", NF90_DOUBLE, [dim_theta], var_theta))
        call nc_check(nf90_def_var(ncid, "zeta", NF90_DOUBLE, [dim_zeta], var_zeta))
        call nc_check(nf90_def_var(ncid, "x", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_x))
        call nc_check(nf90_def_var(ncid, "y", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_y))
        call nc_check(nf90_def_var(ncid, "z", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_z))
        call nc_check(nf90_enddef(ncid))

        call nc_check(nf90_put_var(ncid, var_rho, rho))
        call nc_check(nf90_put_var(ncid, var_theta, theta))
        call nc_check(nf90_put_var(ncid, var_zeta, zeta))
        call nc_check(nf90_put_var(ncid, var_x, x))
        call nc_check(nf90_put_var(ncid, var_y, y))
        call nc_check(nf90_put_var(ncid, var_z, z))
        call nc_check(nf90_close(ncid))
    end subroutine write_missing_units_file

    subroutine write_wrong_dims_file(filename)
        character(len=*), intent(in) :: filename
        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z
        integer, parameter :: nrho = 3, ntheta = 4, nzeta = 5
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp) :: x(nzeta, ntheta, nrho), y(nzeta, ntheta, nrho), z(nzeta, &
                                                                      ntheta, nrho)

        call define_grids(nrho, ntheta, nzeta, rho, theta, zeta)
        x = 0.0_dp
        y = 0.0_dp
        z = 0.0_dp

        call nc_check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "zeta_convention", "cyl"))
        call nc_check(nf90_def_dim(ncid, "rho", nrho, dim_rho))
        call nc_check(nf90_def_dim(ncid, "theta", ntheta, dim_theta))
        call nc_check(nf90_def_dim(ncid, "zeta", nzeta, dim_zeta))
        call nc_check(nf90_def_var(ncid, "rho", NF90_DOUBLE, [dim_rho], var_rho))
        call nc_check(nf90_def_var(ncid, "theta", NF90_DOUBLE, [dim_theta], var_theta))
        call nc_check(nf90_def_var(ncid, "zeta", NF90_DOUBLE, [dim_zeta], var_zeta))
        call nc_check(nf90_def_var(ncid, "x", NF90_DOUBLE, &
                                   [dim_zeta, dim_theta, dim_rho], var_x))
        call nc_check(nf90_def_var(ncid, "y", NF90_DOUBLE, &
                                   [dim_zeta, dim_theta, dim_rho], var_y))
        call nc_check(nf90_def_var(ncid, "z", NF90_DOUBLE, &
                                   [dim_zeta, dim_theta, dim_rho], var_z))
        call nc_check(nf90_put_att(ncid, var_x, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_y, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_z, "units", "cm"))
        call nc_check(nf90_enddef(ncid))

        call nc_check(nf90_put_var(ncid, var_rho, rho))
        call nc_check(nf90_put_var(ncid, var_theta, theta))
        call nc_check(nf90_put_var(ncid, var_zeta, zeta))
        call nc_check(nf90_put_var(ncid, var_x, x))
        call nc_check(nf90_put_var(ncid, var_y, y))
        call nc_check(nf90_put_var(ncid, var_z, z))
        call nc_check(nf90_close(ncid))
    end subroutine write_wrong_dims_file

    subroutine write_zeta_endpoint_file(filename)
        character(len=*), intent(in) :: filename
        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z, &
                   var_num_field_periods
        integer, parameter :: nrho = 3, ntheta = 4, nzeta = 3
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), z(nrho, &
                                                                      ntheta, nzeta)

        call define_grids(nrho, ntheta, nzeta, rho, theta, zeta)
        zeta(nzeta) = TWOPI
        x = 0.0_dp
        y = 0.0_dp
        z = 0.0_dp

        call nc_check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "zeta_convention", "cyl"))
        call nc_check(nf90_def_dim(ncid, "rho", nrho, dim_rho))
        call nc_check(nf90_def_dim(ncid, "theta", ntheta, dim_theta))
        call nc_check(nf90_def_dim(ncid, "zeta", nzeta, dim_zeta))
        call nc_check(nf90_def_var(ncid, "rho", NF90_DOUBLE, [dim_rho], var_rho))
        call nc_check(nf90_def_var(ncid, "theta", NF90_DOUBLE, [dim_theta], var_theta))
        call nc_check(nf90_def_var(ncid, "zeta", NF90_DOUBLE, [dim_zeta], var_zeta))
        call nc_check(nf90_def_var(ncid, "x", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_x))
        call nc_check(nf90_def_var(ncid, "y", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_y))
        call nc_check(nf90_def_var(ncid, "z", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_z))
        call nc_check(nf90_def_var(ncid, "num_field_periods", NF90_INT, &
                                   var_num_field_periods))
        call nc_check(nf90_put_att(ncid, var_x, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_y, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_z, "units", "cm"))
        call nc_check(nf90_enddef(ncid))

        call nc_check(nf90_put_var(ncid, var_rho, rho))
        call nc_check(nf90_put_var(ncid, var_theta, theta))
        call nc_check(nf90_put_var(ncid, var_zeta, zeta))
        call nc_check(nf90_put_var(ncid, var_x, x))
        call nc_check(nf90_put_var(ncid, var_y, y))
        call nc_check(nf90_put_var(ncid, var_z, z))
        call nc_check(nf90_put_var(ncid, var_num_field_periods, 1))
        call nc_check(nf90_close(ncid))
    end subroutine write_zeta_endpoint_file

    subroutine write_unknown_convention_file(filename)
        character(len=*), intent(in) :: filename
        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z, &
                   var_num_field_periods
        integer, parameter :: nrho = 3, ntheta = 4, nzeta = 5
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), z(nrho, &
                                                                      ntheta, nzeta)

        call define_grids(nrho, ntheta, nzeta, rho, theta, zeta)
        x = 0.0_dp
        y = 0.0_dp
        z = 0.0_dp

        call nc_check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "zeta_convention", "unknown"))
        call nc_check(nf90_def_dim(ncid, "rho", nrho, dim_rho))
        call nc_check(nf90_def_dim(ncid, "theta", ntheta, dim_theta))
        call nc_check(nf90_def_dim(ncid, "zeta", nzeta, dim_zeta))
        call nc_check(nf90_def_var(ncid, "rho", NF90_DOUBLE, [dim_rho], var_rho))
        call nc_check(nf90_def_var(ncid, "theta", NF90_DOUBLE, [dim_theta], var_theta))
        call nc_check(nf90_def_var(ncid, "zeta", NF90_DOUBLE, [dim_zeta], var_zeta))
        call nc_check(nf90_def_var(ncid, "x", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_x))
        call nc_check(nf90_def_var(ncid, "y", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_y))
        call nc_check(nf90_def_var(ncid, "z", NF90_DOUBLE, &
                                   [dim_rho, dim_theta, dim_zeta], var_z))
        call nc_check(nf90_def_var(ncid, "num_field_periods", NF90_INT, &
                                   var_num_field_periods))
        call nc_check(nf90_put_att(ncid, var_x, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_y, "units", "cm"))
        call nc_check(nf90_put_att(ncid, var_z, "units", "cm"))
        call nc_check(nf90_enddef(ncid))

        call nc_check(nf90_put_var(ncid, var_rho, rho))
        call nc_check(nf90_put_var(ncid, var_theta, theta))
        call nc_check(nf90_put_var(ncid, var_zeta, zeta))
        call nc_check(nf90_put_var(ncid, var_x, x))
        call nc_check(nf90_put_var(ncid, var_y, y))
        call nc_check(nf90_put_var(ncid, var_z, z))
        call nc_check(nf90_put_var(ncid, var_num_field_periods, 1))
        call nc_check(nf90_close(ncid))
    end subroutine write_unknown_convention_file

end program test_chartmap_validator
