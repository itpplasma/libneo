module chartmap_gate_fixtures
    ! Analytic Boozer-convention chartmap fixtures for the SFL/Boozer chart gate
    ! tests (signed periodic Jacobian, bounded flux-label drift, launch-independent
    ! rotational transform, converged field reconstruction). Self-contained NetCDF
    ! writers in the style of test_chartmap_invert_cart; no external generator.
    !
    ! Circular fixture: concentric circular tori around (CIRC_R_AXIS, 0) with
    ! minor radius r = CIRC_A_MINOR*rho, nfp = 1, and Boozer-like periodic angle
    ! offsets so the chart angles differ from the geometric ones while the rho
    ! surfaces stay exact circles. handedness = +1 makes (rho, theta, zeta)
    ! right-handed (det(dx/du) > 0) by putting theta clockwise in the (R, Z)
    ! plane; handedness = -1 mirrors Z (left-handed chart).
    !
    ! Stellarator fixture: the nfp = 3 rotating-frame construction of
    ! test_chartmap_invert_cart with the same handedness switch.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use math_constants, only: TWOPI
    use netcdf, only: NF90_NETCDF4, NF90_GLOBAL, NF90_DOUBLE, NF90_INT, NF90_NOERR, &
                      nf90_create, nf90_def_dim, nf90_def_var, nf90_put_att, &
                      nf90_enddef, nf90_put_var, nf90_close, nf90_strerror
    implicit none
    private

    public :: write_circular_boozer_chartmap
    public :: write_stellarator_boozer_chartmap
    public :: circular_chart_position

    real(dp), parameter, public :: CIRC_R_AXIS = 180.0_dp
    real(dp), parameter, public :: CIRC_A_MINOR = 45.0_dp
    real(dp), parameter, public :: CIRC_LAM_AMP = 0.04_dp
    real(dp), parameter, public :: CIRC_NU_AMP = 0.03_dp

contains

    subroutine nc_check(status)
        integer, intent(in) :: status
        if (status /= NF90_NOERR) then
            print *, trim(nf90_strerror(status))
            error stop 2
        end if
    end subroutine nc_check

    ! Exact analytic forward map of the circular fixture: chart u = (rho, theta,
    ! zeta) to Cartesian x in cm. Reference truth for the gate tests, independent
    ! of any chart grid resolution.
    subroutine circular_chart_position(rho, theta, zeta, handedness, x)
        real(dp), intent(in) :: rho, theta, zeta
        integer, intent(in) :: handedness
        real(dp), intent(out) :: x(3)

        real(dp) :: r, thg, phg, rmaj

        r = CIRC_A_MINOR*rho
        thg = theta + CIRC_LAM_AMP*rho*sin(theta)
        phg = zeta + CIRC_NU_AMP*rho*sin(thg)
        rmaj = CIRC_R_AXIS + r*cos(thg)
        x(1) = rmaj*cos(phg)
        x(2) = rmaj*sin(phg)
        x(3) = -real(handedness, dp)*r*sin(thg)
    end subroutine circular_chart_position

    subroutine write_circular_boozer_chartmap(path, nrho, ntheta, nzeta, handedness)
        character(len=*), intent(in) :: path
        integer, intent(in) :: nrho, ntheta, nzeta
        integer, intent(in) :: handedness

        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :)
        real(dp) :: pos(3)
        integer :: ir, it, iz

        do ir = 1, nrho
            rho(ir) = real(ir - 1, dp)/real(nrho - 1, dp)
        end do
        do it = 1, ntheta
            theta(it) = TWOPI*real(it - 1, dp)/real(ntheta, dp)
        end do
        do iz = 1, nzeta
            zeta(iz) = TWOPI*real(iz - 1, dp)/real(nzeta, dp)
        end do

        allocate (x(nrho, ntheta, nzeta))
        allocate (y(nrho, ntheta, nzeta))
        allocate (z(nrho, ntheta, nzeta))
        do iz = 1, nzeta
            do it = 1, ntheta
                do ir = 1, nrho
                    call circular_chart_position(rho(ir), theta(it), zeta(iz), &
                                                 handedness, pos)
                    x(ir, it, iz) = pos(1)
                    y(ir, it, iz) = pos(2)
                    z(ir, it, iz) = pos(3)
                end do
            end do
        end do

        call write_chartmap_netcdf(path, nrho, ntheta, nzeta, 1, rho, theta, zeta, &
                                   x, y, z)
    end subroutine write_circular_boozer_chartmap

    ! Reactor-scale synthetic nfp = 3 Boozer chartmap with a genuinely
    ! non-cylindrical toroidal angle (same shaping as test_chartmap_invert_cart).
    subroutine write_stellarator_boozer_chartmap(path, nrho, ntheta, nzeta, handedness)
        character(len=*), intent(in) :: path
        integer, intent(in) :: nrho, ntheta, nzeta
        integer, intent(in) :: handedness

        integer, parameter :: nfp = 3
        real(dp), parameter :: r0 = 180.0_dp
        real(dp), parameter :: a = 45.0_dp
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :)
        real(dp) :: rho_val, theta_val, zeta_val
        real(dp) :: phi_geom, rmaj, zpos, period
        integer :: ir, it, iz

        period = TWOPI/real(nfp, dp)

        do ir = 1, nrho
            rho(ir) = real(ir - 1, dp)/real(nrho - 1, dp)
        end do
        do it = 1, ntheta
            theta(it) = TWOPI*real(it - 1, dp)/real(ntheta, dp)
        end do
        do iz = 1, nzeta
            zeta(iz) = period*real(iz - 1, dp)/real(nzeta, dp)
        end do

        allocate (x(nrho, ntheta, nzeta))
        allocate (y(nrho, ntheta, nzeta))
        allocate (z(nrho, ntheta, nzeta))
        do iz = 1, nzeta
            zeta_val = zeta(iz)
            do it = 1, ntheta
                theta_val = theta(it)
                do ir = 1, nrho
                    rho_val = rho(ir)
                    phi_geom = zeta_val + 0.08_dp*rho_val*sin(theta_val) + &
                               0.02_dp*sin(2.0_dp*zeta_val)
                    rmaj = r0 + a*rho_val*cos(theta_val) + &
                           0.04_dp*a*rho_val*cos(theta_val - 2.0_dp*zeta_val)
                    zpos = a*rho_val*sin(theta_val) + &
                           0.02_dp*a*rho_val*sin(zeta_val)

                    x(ir, it, iz) = rmaj*cos(phi_geom)
                    y(ir, it, iz) = rmaj*sin(phi_geom)
                    z(ir, it, iz) = -real(handedness, dp)*zpos
                end do
            end do
        end do

        call write_chartmap_netcdf(path, nrho, ntheta, nzeta, nfp, rho, theta, zeta, &
                                   x, y, z)
    end subroutine write_stellarator_boozer_chartmap

    subroutine write_chartmap_netcdf(path, nrho, ntheta, nzeta, nfp, rho, theta, &
                                     zeta, x, y, z)
        character(len=*), intent(in) :: path
        integer, intent(in) :: nrho, ntheta, nzeta, nfp
        real(dp), intent(in) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp), intent(in) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), &
                                z(nrho, ntheta, nzeta)

        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z, var_nfp

        call nc_check(nf90_create(trim(path), NF90_NETCDF4, ncid))
        call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "zeta_convention", "boozer"))
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
        call nc_check(nf90_def_var(ncid, "num_field_periods", NF90_INT, var_nfp))
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
        call nc_check(nf90_put_var(ncid, var_nfp, nfp))
        call nc_check(nf90_close(ncid))
    end subroutine write_chartmap_netcdf

end module chartmap_gate_fixtures
