program test_refcoords_file_detection
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: detect_refcoords_file_type, &
                                  refcoords_file_chartmap, refcoords_file_vmec_wout
    use math_constants, only: TWOPI
    use netcdf
    implicit none

    integer :: ierr
    integer :: file_type
    character(len=2048) :: message
    integer :: nerrors

    nerrors = 0

    call detect_refcoords_file_type("chartmap.nc", file_type, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: detect_refcoords_file_type(chartmap.nc): ", trim(message)
        nerrors = nerrors + 1
    else if (file_type /= refcoords_file_chartmap) then
        print *, "  FAIL: expected CHARTMAP file_type for chartmap.nc, got ", file_type
        nerrors = nerrors + 1
    else
        print *, "  PASS: chartmap.nc detected as CHARTMAP"
    end if

    call detect_refcoords_file_type("wout.nc", file_type, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: detect_refcoords_file_type(wout.nc): ", trim(message)
        nerrors = nerrors + 1
    else if (file_type /= refcoords_file_vmec_wout) then
        print *, "  FAIL: expected VMEC_WOUT file_type for wout.nc, got ", file_type
        nerrors = nerrors + 1
    else
        print *, "  PASS: wout.nc detected as VMEC_WOUT"
    end if

    call write_invalid_chartmap_missing_convention("chartmap_missing_convention.nc")
    call detect_refcoords_file_type("chartmap_missing_convention.nc", file_type, &
                                    ierr, message)
    if (file_type /= refcoords_file_chartmap) then
        print *, &
            "  FAIL: expected CHARTMAP file_type for missing convention file, got ", &
            file_type
        nerrors = nerrors + 1
    else if (ierr == 0 .or. index(message, "zeta_convention") == 0) then
        print *, "  FAIL: expected zeta_convention error for missing convention file"
        print *, trim(message)
        nerrors = nerrors + 1
    else
        print *, "  PASS: invalid chartmap detected and reports zeta_convention error"
    end if

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) detected in refcoords file detection"
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

    subroutine write_invalid_chartmap_missing_convention(filename)
        character(len=*), intent(in) :: filename

        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta
        integer :: var_x, var_y, var_z
        integer :: var_num_field_periods
        integer, parameter :: nrho = 3, ntheta = 4, nzeta = 5
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), z(nrho, ntheta, &
                                                                      nzeta)
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

        x = 0.0_dp
        y = 0.0_dp
        z = 0.0_dp

        call nc_check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
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
    end subroutine write_invalid_chartmap_missing_convention

end program test_refcoords_file_detection
