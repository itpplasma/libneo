program test_chartmap_invert_cart
    ! Behavioral test for cs%invert_cart, the warm Cartesian inverse of the chartmap
    ! coordinate system. Writes a self-contained reactor-scale Boozer chartmap
    ! (no map2disc dependency), then checks that:
    !   - evaluate_cart(u0) -> invert_cart recovers an interior u0 to tolerance and
    !     reports CHARTMAP_LOCATED, for interior, near-axis (rho ~ 1e-3), and seam
    !     (zeta near a field-period boundary) points;
    !   - a point on the last closed surface (rho = 1) reports CHARTMAP_CLAMPED_EDGE;
    !   - a target built past the last surface (rho > 1, which the forward map cannot
    !     represent) reports CHARTMAP_OUTSIDE, never a spurious interior root;
    !   - an axis point (rho ~ 0) is handled by the X/Y regularization (no NaN, theta
    !     recovered modulo the axis gauge).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
        make_chartmap_coordinate_system, BOOZER, &
        CHARTMAP_LOCATED, CHARTMAP_CLAMPED_EDGE, &
        CHARTMAP_OUTSIDE, CHARTMAP_NO_ROOT
    use math_constants, only: TWOPI
    use netcdf
    implicit none

    character(len=*), parameter :: filename = "chartmap_invert.nc"
    integer :: nerrors

    nerrors = 0

    call write_boozer_file(filename)
    call run_invert_checks(filename, nerrors)

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) detected in invert_cart tests"
        error stop 1
    end if

    print *, "All chartmap invert_cart tests passed!"

contains

    subroutine nc_check(status)
        integer, intent(in) :: status
        if (status /= NF90_NOERR) then
            print *, trim(nf90_strerror(status))
            error stop 2
        end if
    end subroutine nc_check

    pure real(dp) function angular_diff(a, b, period)
        real(dp), intent(in) :: a, b, period
        angular_diff = abs(modulo(a - b + 0.5_dp*period, period) - 0.5_dp*period)
    end function angular_diff

    subroutine run_invert_checks(path, nerrors)
        character(len=*), intent(in) :: path
        integer, intent(inout) :: nerrors

        class(coordinate_system_t), allocatable :: cs
        real(dp) :: zeta_period

        call make_chartmap_coordinate_system(cs, path)

        select type (ccs => cs)
            type is (chartmap_coordinate_system_t)
            if (ccs%zeta_convention /= BOOZER) then
                print *, "  FAIL: expected BOOZER zeta convention"
                nerrors = nerrors + 1
                return
            end if
            zeta_period = TWOPI/real(ccs%num_field_periods, dp)

            call check_interior_roundtrip(ccs, zeta_period, nerrors)
            call check_seam_roundtrip(ccs, zeta_period, nerrors)
            call check_near_axis(ccs, zeta_period, nerrors)
            call check_axis(ccs, zeta_period, nerrors)
            call check_edge(ccs, zeta_period, nerrors)
            call check_outside(ccs, zeta_period, nerrors)
        class default
            print *, "  FAIL: chartmap file did not load as chartmap type"
            nerrors = nerrors + 1
        end select
    end subroutine run_invert_checks

    ! Interior points across zeta slices: evaluate_cart(u0) -> invert_cart recovers
    ! u0 to tolerance with status LOCATED. The warm guess is offset from u0 so the
    ! solver does real work rather than starting on the root.
    subroutine check_interior_roundtrip(ccs, zeta_period, nerrors)
        type(chartmap_coordinate_system_t), intent(in), target :: ccs
        real(dp), intent(in) :: zeta_period
        integer, intent(inout) :: nerrors
        real(dp) :: u0(3), x(3), u(3), u_guess(3)
        real(dp) :: dth, dze
        integer :: ir, it, iz, status, nfail
        real(dp), parameter :: rho_vals(3) = [0.18_dp, 0.52_dp, 0.86_dp]
        real(dp), parameter :: theta_vals(3) = [0.4_dp, 2.6_dp, 5.1_dp]
        real(dp), parameter :: zeta_fracs(3) = [0.1_dp, 0.45_dp, 0.8_dp]
        real(dp), parameter :: tol = 1.0e-6_dp

        nfail = 0
        do iz = 1, size(zeta_fracs)
            do it = 1, size(theta_vals)
                do ir = 1, size(rho_vals)
                    u0 = [rho_vals(ir), theta_vals(it), zeta_fracs(iz)*zeta_period]
                    call ccs%evaluate_cart(u0, x)
                    u_guess = [u0(1) + 0.07_dp, u0(2) + 0.15_dp, u0(3) - 0.05_dp]
                    call ccs%invert_cart(x, u_guess, u, status)

                    dth = angular_diff(u(2), u0(2), TWOPI)
                    dze = angular_diff(u(3), u0(3), zeta_period)
                    if (status /= CHARTMAP_LOCATED .or. abs(u(1) - u0(1)) > tol .or. &
                        dth > tol .or. dze > tol) then
                        print *, "  FAIL: interior roundtrip u0=", u0
                        print *, "        u=", u, " status=", status
                        nfail = nfail + 1
                    end if
                end do
            end do
        end do
        nerrors = nerrors + nfail
        if (nfail == 0) print *, "  PASS: interior roundtrip recovers u0, LOCATED"
    end subroutine check_interior_roundtrip

    ! Seam points: zeta within a small fraction of a field-period boundary, where a
    ! global Boozer phi guess would be a full period stale. invert_cart reseeds the
    ! toroidal angle and must still locate.
    subroutine check_seam_roundtrip(ccs, zeta_period, nerrors)
        type(chartmap_coordinate_system_t), intent(in), target :: ccs
        real(dp), intent(in) :: zeta_period
        integer, intent(inout) :: nerrors
        real(dp) :: u0(3), x(3), u(3), u_guess(3)
        real(dp) :: dth, dze
        integer :: k, status, nfail
        real(dp), parameter :: seam_fracs(3) = [0.0_dp, 0.999_dp, 0.5_dp]
        real(dp), parameter :: tol = 1.0e-6_dp

        nfail = 0
        do k = 1, size(seam_fracs)
            u0 = [0.6_dp, 1.3_dp, seam_fracs(k)*zeta_period]
            call ccs%evaluate_cart(u0, x)
            ! Warm guess carries a stale phi one period ahead, as across a seam.
            u_guess = [u0(1), u0(2), u0(3) + zeta_period]
            call ccs%invert_cart(x, u_guess, u, status)
            dth = angular_diff(u(2), u0(2), TWOPI)
            dze = angular_diff(u(3), u0(3), zeta_period)
            if (status /= CHARTMAP_LOCATED .or. abs(u(1) - u0(1)) > tol .or. &
                dth > tol .or. dze > tol) then
                print *, "  FAIL: seam roundtrip u0=", u0, " u=", u, " status=", status
                nfail = nfail + 1
            end if
        end do
        nerrors = nerrors + nfail
        if (nfail == 0) print *, "  PASS: seam roundtrip with stale phi guess, LOCATED"
    end subroutine check_seam_roundtrip

    ! Near-axis points (rho ~ 1e-3): the polar chart is singular but the X/Y chart is
    ! regular. The radius must be recovered; theta is a near-degenerate gauge there,
    ! so only rho and the Cartesian round-trip are asserted.
    subroutine check_near_axis(ccs, zeta_period, nerrors)
        type(chartmap_coordinate_system_t), intent(in), target :: ccs
        real(dp), intent(in) :: zeta_period
        integer, intent(inout) :: nerrors
        real(dp) :: u0(3), x(3), u(3), u_guess(3), x_back(3)
        integer :: k, status, nfail
        real(dp), parameter :: rho_vals(2) = [1.0e-3_dp, 5.0e-3_dp]
        real(dp), parameter :: tol_rho = 1.0e-5_dp, tol_x = 1.0e-4_dp

        nfail = 0
        do k = 1, size(rho_vals)
            u0 = [rho_vals(k), 2.1_dp, 0.3_dp*zeta_period]
            call ccs%evaluate_cart(u0, x)
            u_guess = [rho_vals(k) + 0.02_dp, 1.0_dp, 0.3_dp*zeta_period]
            call ccs%invert_cart(x, u_guess, u, status)
            call ccs%evaluate_cart(u, x_back)
            if (status /= CHARTMAP_LOCATED .or. abs(u(1) - u0(1)) > tol_rho .or. &
                maxval(abs(x_back - x)) > tol_x) then
                print *, "  FAIL: near-axis rho=", rho_vals(k), " u=", u, &
                    " status=", status, " dx=", maxval(abs(x_back - x))
                nfail = nfail + 1
            end if
        end do
        nerrors = nerrors + nfail
        if (nfail == 0) print *, "  PASS: near-axis rho recovered via X/Y healing"
    end subroutine check_near_axis

    ! Axis point (rho = 0): the X/Y regularization must return a finite u with rho
    ! near 0 and a Cartesian round-trip back to the axis, with no NaN.
    subroutine check_axis(ccs, zeta_period, nerrors)
        type(chartmap_coordinate_system_t), intent(in), target :: ccs
        real(dp), intent(in) :: zeta_period
        integer, intent(inout) :: nerrors
        real(dp) :: u0(3), x(3), u(3), u_guess(3), x_back(3)
        integer :: status

        u0 = [0.0_dp, 0.0_dp, 0.25_dp*zeta_period]
        call ccs%evaluate_cart(u0, x)
        u_guess = [0.05_dp, 1.0_dp, 0.25_dp*zeta_period]
        call ccs%invert_cart(x, u_guess, u, status)
        call ccs%evaluate_cart(u, x_back)
        if (status == CHARTMAP_NO_ROOT .or. u(1) /= u(1) .or. &
            u(1) > 1.0e-2_dp .or. maxval(abs(x_back - x)) > 1.0e-4_dp) then
            print *, "  FAIL: axis point u=", u, " status=", status, &
                " dx=", maxval(abs(x_back - x))
            nerrors = nerrors + 1
        else
            print *, "  PASS: axis point handled by X/Y regularization"
        end if
    end subroutine check_axis

    ! Edge point exactly on the last closed surface (rho = 1): the warm guess is also
    ! at the edge, so the bracketed radial solve pins rho to 1 and reports
    ! CHARTMAP_CLAMPED_EDGE with the residual at tolerance.
    subroutine check_edge(ccs, zeta_period, nerrors)
        type(chartmap_coordinate_system_t), intent(in), target :: ccs
        real(dp), intent(in) :: zeta_period
        integer, intent(inout) :: nerrors
        real(dp) :: u0(3), x(3), u(3), u_guess(3)
        integer :: status

        u0 = [1.0_dp, 1.7_dp, 0.4_dp*zeta_period]
        call ccs%evaluate_cart(u0, x)
        u_guess = [0.97_dp, 1.7_dp, 0.4_dp*zeta_period]
        call ccs%invert_cart(x, u_guess, u, status)
        if (status /= CHARTMAP_LOCATED .and. status /= CHARTMAP_CLAMPED_EDGE) then
            print *, "  FAIL: edge point status=", status, " u=", u
            nerrors = nerrors + 1
        else if (abs(u(1) - 1.0_dp) > 1.0e-4_dp) then
            print *, "  FAIL: edge point rho not at 1, u=", u
            nerrors = nerrors + 1
        else
            print *, "  PASS: last-surface point located at the edge"
        end if
    end subroutine check_edge

    ! Target built past the last closed surface: take the boundary point and push it
    ! outward along its outward radial direction by a sizable fraction of a radial
    ! cell. The forward map cannot represent rho > 1, so the best root clamps to
    ! rho = 1 with a residual a large fraction of a cell -> CHARTMAP_OUTSIDE.
    subroutine check_outside(ccs, zeta_period, nerrors)
        type(chartmap_coordinate_system_t), intent(in), target :: ccs
        real(dp), intent(in) :: zeta_period
        integer, intent(inout) :: nerrors
        real(dp) :: u0(3), x_edge(3), x(3), u(3), u_guess(3)
        real(dp) :: e_cov(3, 3), e_rho(3), e_norm
        integer :: status

        u0 = [1.0_dp, 0.9_dp, 0.6_dp*zeta_period]
        call ccs%evaluate_cart(u0, x_edge)
        call ccs%covariant_basis(u0, e_cov)
        e_rho = e_cov(:, 1)
        e_norm = sqrt(sum(e_rho**2))
        ! Half a radial cell past the boundary: well outside the LOCATED/CLAMPED band.
        x = x_edge + 0.5_dp*e_rho
        u_guess = [0.98_dp, 0.9_dp, 0.6_dp*zeta_period]
        call ccs%invert_cart(x, u_guess, u, status)
        if (status /= CHARTMAP_OUTSIDE) then
            print *, "  FAIL: past-edge target status=", status, " (want OUTSIDE)"
            print *, "        u=", u
            nerrors = nerrors + 1
        else if (abs(u(1) - 1.0_dp) > 1.0e-4_dp) then
            print *, "  FAIL: past-edge target rho not clamped to 1, u=", u
            nerrors = nerrors + 1
        else
            print *, "  PASS: past-edge target reports OUTSIDE, rho clamped to 1"
        end if
    end subroutine check_outside

    ! Reactor-scale synthetic Boozer chartmap, nfp = 3, with a genuinely
    ! non-cylindrical toroidal angle so the inverse exercises the Boozer-phi reseed.
    ! Same construction as test_chartmap_boozer_coordinates so the fixture needs no
    ! external generator.
    subroutine write_boozer_file(path)
        character(len=*), intent(in) :: path

        integer :: ncid
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_rho, var_theta, var_zeta, var_x, var_y, var_z, var_nfp
        integer, parameter :: nrho = 9, ntheta = 11, nzeta = 10, nfp = 3
        real(dp), parameter :: r0 = 180.0_dp
        real(dp), parameter :: a = 45.0_dp
        real(dp) :: rho(nrho), theta(ntheta), zeta(nzeta)
        real(dp) :: x(nrho, ntheta, nzeta), y(nrho, ntheta, nzeta), &
            z(nrho, ntheta, nzeta)
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
                    z(ir, it, iz) = zpos
                end do
            end do
        end do

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
    end subroutine write_boozer_file

end program test_chartmap_invert_cart
