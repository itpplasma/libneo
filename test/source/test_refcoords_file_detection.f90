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

    call write_good_chartmap_file("chartmap_refcoords.nc")
    call detect_refcoords_file_type("chartmap_refcoords.nc", file_type, ierr, &
                                    message)
    if (ierr /= 0) then
        print *, "  FAIL: detect_refcoords_file_type(chartmap_refcoords.nc): ", &
            trim(message)
        nerrors = nerrors + 1
    else if (file_type /= refcoords_file_chartmap) then
        print *, &
            "  FAIL: expected CHARTMAP file_type for chartmap_refcoords.nc, got ", &
            file_type
        nerrors = nerrors + 1
    else
        print *, "  PASS: chartmap_refcoords.nc detected as CHARTMAP"
    end if

    call write_good_boozer_file("chartmap_refcoords_boozer.nc")
    call detect_refcoords_file_type("chartmap_refcoords_boozer.nc", file_type, ierr, &
                                    message)
    if (ierr /= 0) then
        print *, "  FAIL: detect_refcoords_file_type(chartmap_refcoords_boozer.nc): ", &
            trim(message)
        nerrors = nerrors + 1
    else if (file_type /= refcoords_file_chartmap) then
        print *, &
            "  FAIL: expected CHARTMAP file_type for chartmap_refcoords_boozer.nc, got ", &
            file_type
        nerrors = nerrors + 1
    else
        print *, "  PASS: chartmap_refcoords_boozer.nc detected as CHARTMAP"
    end if

    call write_good_wout_file("wout_refcoords.nc")
    call detect_refcoords_file_type("wout_refcoords.nc", file_type, ierr, message)
    if (ierr /= 0) then
        print *, "  FAIL: detect_refcoords_file_type(wout_refcoords.nc): ", &
            trim(message)
        nerrors = nerrors + 1
    else if (file_type /= refcoords_file_vmec_wout) then
        print *, &
            "  FAIL: expected VMEC_WOUT file_type for wout_refcoords.nc, got ", &
            file_type
        nerrors = nerrors + 1
    else
        print *, "  PASS: wout_refcoords.nc detected as VMEC_WOUT"
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

    subroutine write_good_chartmap_file(filename)
        character(len=*), intent(in) :: filename

        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z
        integer :: var_num_field_periods
        integer, parameter :: nrho = 3, ntheta = 4, nzeta = 5
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), z(nrho, &
                                                                      ntheta, nzeta)
        real(dp) :: rho_val, theta_val, zeta_val, rmaj, zpos, period
        integer :: i, ir, it, iz

        period = TWOPI

        do i = 1, nrho
            rho(i) = real(i - 1, dp)/real(nrho - 1, dp)
        end do
        do i = 1, ntheta
            theta(i) = TWOPI*real(i - 1, dp)/real(ntheta, dp)
        end do
        do i = 1, nzeta
            zeta(i) = period*real(i - 1, dp)/real(nzeta, dp)
        end do

        do iz = 1, nzeta
            zeta_val = zeta(iz)
            do it = 1, ntheta
                theta_val = theta(it)
                do ir = 1, nrho
                    rho_val = rho(ir)
                    rmaj = 150.0_dp + 35.0_dp*rho_val*cos(theta_val)
                    zpos = 35.0_dp*rho_val*sin(theta_val)

                    x(ir, it, iz) = rmaj*cos(zeta_val)
                    y(ir, it, iz) = rmaj*sin(zeta_val)
                    z(ir, it, iz) = zpos
                end do
            end do
        end do

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
    end subroutine write_good_chartmap_file

    subroutine write_good_wout_file(filename)
        character(len=*), intent(in) :: filename

        integer :: ncid
        integer :: var_rmnc

        call nc_check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
        call nc_check(nf90_def_var(ncid, "rmnc", NF90_DOUBLE, var_rmnc))
        call nc_check(nf90_put_var(ncid, var_rmnc, 0.0_dp))
        call nc_check(nf90_close(ncid))
    end subroutine write_good_wout_file

    subroutine write_good_boozer_file(filename)
        character(len=*), intent(in) :: filename

        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z
        integer :: var_num_field_periods
        integer, parameter :: nrho = 3, ntheta = 5, nzeta = 6, nfp = 2
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), z(nrho, &
                                                                      ntheta, nzeta)
        real(dp) :: rho_val, theta_val, zeta_val, phi_geom, rmaj, zpos, period
        integer :: i, ir, it, iz

        period = TWOPI/real(nfp, dp)

        do i = 1, nrho
            rho(i) = real(i - 1, dp)/real(nrho - 1, dp)
        end do
        do i = 1, ntheta
            theta(i) = TWOPI*real(i - 1, dp)/real(ntheta, dp)
        end do
        do i = 1, nzeta
            zeta(i) = period*real(i - 1, dp)/real(nzeta, dp)
        end do

        do iz = 1, nzeta
            zeta_val = zeta(iz)
            do it = 1, ntheta
                theta_val = theta(it)
                do ir = 1, nrho
                    rho_val = rho(ir)
                    phi_geom = zeta_val + 0.08_dp*rho_val*sin(theta_val)
                    rmaj = 150.0_dp + 35.0_dp*rho_val*cos(theta_val)
                    zpos = 35.0_dp*rho_val*sin(theta_val)

                    x(ir, it, iz) = rmaj*cos(phi_geom)
                    y(ir, it, iz) = rmaj*sin(phi_geom)
                    z(ir, it, iz) = zpos
                end do
            end do
        end do

        call nc_check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "zeta_convention", "boozer"))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "rho_convention", "rho_tor"))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "boozer_field", 1))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "rmajor", 150.0_dp))
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
        call nc_check(nf90_put_var(ncid, var_num_field_periods, nfp))
        call nc_check(nf90_close(ncid))
    end subroutine write_good_boozer_file

end program test_refcoords_file_detection
