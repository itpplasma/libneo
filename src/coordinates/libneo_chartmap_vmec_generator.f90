module libneo_chartmap_vmec_generator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: vmec_coordinate_system_t
    use math_constants, only: TWOPI
    use netcdf
    use new_vmec_stuff_mod, only: netcdffile, nper
    use spline_vmec_sub, only: spline_vmec_data
    implicit none

contains

    subroutine write_chartmap_from_vmec(wout_file, outfile, nrho, ntheta, nzeta, &
                                        ierr, message)
        character(len=*), intent(in) :: wout_file
        character(len=*), intent(in) :: outfile
        integer, intent(in) :: nrho
        integer, intent(in) :: ntheta
        integer, intent(in) :: nzeta
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        real(dp), allocatable :: rho(:), theta(:), zeta(:)
        real(dp), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :)

        ierr = 0
        message = ""

        if (nrho < 2 .or. ntheta < 2 .or. nzeta < 1) then
            ierr = 1
            message = "invalid grid size(s)"
            return
        end if
        if (nrho > 65 .or. ntheta > 65 .or. nzeta > 65) then
            ierr = 1
            message = "grid size(s) must be <= 65"
            return
        end if

        netcdffile = trim(wout_file)
        call spline_vmec_data
        if (nper < 1) then
            ierr = 2
            message = "VMEC nper invalid"
            return
        end if

        allocate (rho(nrho), theta(ntheta), zeta(nzeta))
        call build_grids(rho, theta, zeta, nper)

        allocate (x(nrho, ntheta, nzeta))
        allocate (y(nrho, ntheta, nzeta))
        allocate (z(nrho, ntheta, nzeta))
        call sample_vmec_cartesian(rho, theta, zeta, x, y, z)

        call write_chartmap_file(outfile, rho, theta, zeta, x, y, z, nper, &
                                 ierr, message)
    end subroutine write_chartmap_from_vmec

    subroutine build_grids(rho, theta, zeta, num_field_periods)
        real(dp), intent(out) :: rho(:)
        real(dp), intent(out) :: theta(:)
        real(dp), intent(out) :: zeta(:)
        integer, intent(in) :: num_field_periods

        integer :: i
        real(dp) :: zeta_period

        zeta_period = TWOPI/real(num_field_periods, dp)

        do i = 1, size(rho)
            rho(i) = real(i - 1, dp)/real(size(rho) - 1, dp)
        end do
        do i = 1, size(theta)
            theta(i) = TWOPI*real(i - 1, dp)/real(size(theta), dp)
        end do
        do i = 1, size(zeta)
            zeta(i) = zeta_period*real(i - 1, dp)/real(size(zeta), dp)
        end do
    end subroutine build_grids

    subroutine sample_vmec_cartesian(rho, theta, zeta, x, y, z)
        real(dp), intent(in) :: rho(:)
        real(dp), intent(in) :: theta(:)
        real(dp), intent(in) :: zeta(:)
        real(dp), intent(out) :: x(:, :, :)
        real(dp), intent(out) :: y(:, :, :)
        real(dp), intent(out) :: z(:, :, :)

        type(vmec_coordinate_system_t) :: vmec
        real(dp) :: u(3), xcyl(3)
        real(dp) :: R, ph, Zc
        integer :: ir, it, iz

        do iz = 1, size(zeta)
            ph = zeta(iz)
            do it = 1, size(theta)
                u(2) = theta(it)
                do ir = 1, size(rho)
                    u(1) = rho(ir)**2
                    u(3) = ph
                    call vmec%evaluate_cyl(u, xcyl)
                    R = xcyl(1)
                    Zc = xcyl(3)
                    x(ir, it, iz) = R*cos(ph)
                    y(ir, it, iz) = R*sin(ph)
                    z(ir, it, iz) = Zc
                end do
            end do
        end do
    end subroutine sample_vmec_cartesian

    subroutine write_chartmap_file(outfile, rho, theta, zeta, x, y, z, &
                                   num_field_periods, &
                                   ierr, message)
        character(len=*), intent(in) :: outfile
        real(dp), intent(in) :: rho(:)
        real(dp), intent(in) :: theta(:)
        real(dp), intent(in) :: zeta(:)
        real(dp), intent(in) :: x(:, :, :)
        real(dp), intent(in) :: y(:, :, :)
        real(dp), intent(in) :: z(:, :, :)
        integer, intent(in) :: num_field_periods
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta
        integer :: var_x, var_y, var_z, var_num_field_periods

        ierr = 0
        message = ""

        call nc_check(nf90_create(trim(outfile), NF90_NETCDF4, ncid), ierr, message)
        if (ierr /= 0) return

        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "zeta_convention", "cyl"), &
                      ierr, message)
        if (ierr /= 0) return

        call nc_check(nf90_def_dim(ncid, "rho", size(rho), dim_rho), ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_def_dim(ncid, "theta", size(theta), dim_theta), &
                      ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_def_dim(ncid, "zeta", size(zeta), dim_zeta), ierr, message)
        if (ierr /= 0) return

        call nc_check(nf90_def_var(ncid, "rho", NF90_DOUBLE, [dim_rho], var_rho), &
                      ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_def_var(ncid, "theta", NF90_DOUBLE, [dim_theta], &
                                   var_theta), &
                      ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_def_var(ncid, "zeta", NF90_DOUBLE, [dim_zeta], var_zeta), &
                      ierr, message)
        if (ierr /= 0) return

        call nc_check(nf90_def_var(ncid, "x", NF90_DOUBLE, [dim_rho, dim_theta, &
                                                            dim_zeta], &
                                   var_x), ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_def_var(ncid, "y", NF90_DOUBLE, [dim_rho, dim_theta, &
                                                            dim_zeta], &
                                   var_y), ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_def_var(ncid, "z", NF90_DOUBLE, [dim_rho, dim_theta, &
                                                            dim_zeta], &
                                   var_z), ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_def_var(ncid, "num_field_periods", NF90_INT, &
                                   var_num_field_periods), ierr, message)
        if (ierr /= 0) return

        call nc_check(nf90_put_att(ncid, var_x, "units", "cm"), ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_put_att(ncid, var_y, "units", "cm"), ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_put_att(ncid, var_z, "units", "cm"), ierr, message)
        if (ierr /= 0) return

        call nc_check(nf90_enddef(ncid), ierr, message)
        if (ierr /= 0) return

        call nc_check(nf90_put_var(ncid, var_rho, rho), ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_put_var(ncid, var_theta, theta), ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_put_var(ncid, var_zeta, zeta), ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_put_var(ncid, var_x, x), ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_put_var(ncid, var_y, y), ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_put_var(ncid, var_z, z), ierr, message)
        if (ierr /= 0) return
        call nc_check(nf90_put_var(ncid, var_num_field_periods, num_field_periods), &
                      ierr, message)
        if (ierr /= 0) return

        call nc_check(nf90_close(ncid), ierr, message)
    end subroutine write_chartmap_file

    subroutine nc_check(status, ierr, message)
        integer, intent(in) :: status
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        ierr = 0
        message = ""
        if (status /= NF90_NOERR) then
            ierr = 1
            message = trim(nf90_strerror(status))
        end if
    end subroutine nc_check

end module libneo_chartmap_vmec_generator
