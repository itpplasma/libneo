program test_chartmap_gate_fieldline
    ! Boozer/SFL chart gates (b) and (c) (ROADMAP straight-field-line gauge gates):
    !   (b) bounded flux-label drift over N toroidal transits of a traced field
    !       line pulled back through the chart inverse;
    !   (c) launch-independent rotational transform from an ensemble of launch
    !       points and angles covering radii, poloidal angles, and toroidal
    !       sectors of the full torus.
    !
    ! Field: the production circular tokamak field (neo_circular_tokamak_field),
    ! B_R = -Bp (Z - Z0)/R, B_phi = Bt, B_Z = Bp (R - R0)/R. Its field lines lie
    ! exactly on circular tori r = const around (R0, Z0) and are straight in the
    ! geometric angles: alpha(phi) = alpha0 + (Bp/Bt) (phi - phi0). The trace
    ! therefore uses the closed-form line (zero integrator noise); a tangency
    ! self-check first verifies at machine precision that the closed form is a
    ! field line of the production evaluator.
    !
    ! Chart: the analytic circular Boozer-convention fixture whose rho surfaces
    ! are exactly the field's flux surfaces (r = CIRC_A_MINOR*rho), with periodic
    ! Boozer-like angle offsets so chart angles differ from geometric ones. The
    ! fixture is right-handed (theta clockwise in the (R, Z) plane), so the
    ! transform measured in chart angles is iota_chart = -Bp/Bt.
    !
    ! Every pullback uses the warm single-point inverse invert_cart with the
    ! previous sample as guess -- the orbit-pusher usage pattern.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
                                  make_chartmap_coordinate_system, CHARTMAP_LOCATED
    use math_constants, only: TWOPI
    use neo_circular_tokamak_field, only: circular_tokamak_field_t
    use chartmap_gate_fixtures, only: write_circular_boozer_chartmap, &
                                      circular_chart_position, CIRC_R_AXIS, &
                                      CIRC_A_MINOR
    implicit none

    character(len=*), parameter :: chart_path = "chartmap_gate_fieldline.nc"
    real(dp), parameter :: B_POL = 0.3183099_dp
    real(dp), parameter :: B_TOR = 1.0_dp
    real(dp), parameter :: iota_geom = B_POL/B_TOR
    real(dp), parameter :: iota_chart_ref = -iota_geom
    integer, parameter :: n_transit = 50
    integer, parameter :: n_per_transit = 12
    real(dp), parameter :: tol_drift = 1.0e-5_dp
    real(dp), parameter :: tol_iota_abs = 5.0e-4_dp
    real(dp), parameter :: tol_iota_spread = 5.0e-4_dp

    type(circular_tokamak_field_t) :: field
    class(coordinate_system_t), allocatable :: cs
    integer :: nerrors

    nerrors = 0

    call field%circular_tokamak_field_init(R_axis=CIRC_R_AXIS, Z_axis=0.0_dp, &
                                           B_pol_ampl=B_POL, B_tor_ampl=B_TOR)
    call check_closed_form_is_field_line(field, nerrors)

    call write_circular_boozer_chartmap(chart_path, 17, 48, 12, 1)
    call make_chartmap_coordinate_system(cs, chart_path)

    select type (ccs => cs)
    type is (chartmap_coordinate_system_t)
        call run_ensemble_gates(ccs, nerrors)
    class default
        print *, "  FAIL: fixture did not load as chartmap type"
        nerrors = nerrors + 1
    end select

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) detected in field-line gate tests"
        error stop 1
    end if

    print *, "All chartmap field-line gate tests passed!"

contains

    ! Verify the closed-form line alpha(phi) = alpha0 + iota_geom*(phi - phi0) is
    ! a field line of the production evaluator: the Cartesian curve tangent must
    ! be parallel to B (relative mismatch at rounding level).
    subroutine check_closed_form_is_field_line(fld, nerrors)
        type(circular_tokamak_field_t), intent(in) :: fld
        integer, intent(inout) :: nerrors

        real(dp), parameter :: tol = 1.0e-12_dp
        real(dp) :: r, alpha, phi, rmaj, zpos
        real(dp) :: dr_dphi, dz_dphi, tangent(3), b_cyl(3), b_cart(3), scale
        real(dp) :: rel
        integer :: k, m, nfail

        nfail = 0
        do m = 1, 3
            r = CIRC_A_MINOR*(0.2_dp + 0.3_dp*real(m - 1, dp))
            do k = 1, 7
                alpha = 0.83_dp*real(k, dp)
                phi = 1.7_dp*real(k, dp)
                rmaj = CIRC_R_AXIS + r*cos(alpha)
                zpos = r*sin(alpha)
                dr_dphi = -iota_geom*r*sin(alpha)
                dz_dphi = iota_geom*r*cos(alpha)
                tangent = [dr_dphi*cos(phi) - rmaj*sin(phi), &
                           dr_dphi*sin(phi) + rmaj*cos(phi), dz_dphi]
                call fld%compute_bfield([rmaj, phi, zpos], b_cyl)
                b_cart = [b_cyl(1)*cos(phi) - b_cyl(2)*sin(phi), &
                          b_cyl(1)*sin(phi) + b_cyl(2)*cos(phi), b_cyl(3)]
                scale = rmaj/b_cyl(2)
                rel = sqrt(sum((tangent - scale*b_cart)**2))/sqrt(sum(tangent**2))
                if (rel > tol) then
                    print *, "  FAIL: closed form not tangent to B, rel=", rel
                    nfail = nfail + 1
                end if
            end do
        end do
        nerrors = nerrors + nfail
        if (nfail == 0) print *, &
            "  PASS: closed-form trace is a field line of the production evaluator"
    end subroutine check_closed_form_is_field_line

    subroutine run_ensemble_gates(ccs, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        integer, intent(inout) :: nerrors

        real(dp), parameter :: rho_vals(3) = [0.30_dp, 0.55_dp, 0.80_dp]
        real(dp), parameter :: theta_vals(3) = [0.5_dp, 2.6_dp, 4.9_dp]
        real(dp), parameter :: zeta_vals(3) = [0.3_dp, 2.2_dp, 4.5_dp]
        real(dp) :: iota_est(27), drift(27)
        real(dp) :: iota_min, iota_max
        integer :: ir, it, iz, m, nfail_trace, nfail

        m = 0
        nfail_trace = 0
        do iz = 1, size(zeta_vals)
            do it = 1, size(theta_vals)
                do ir = 1, size(rho_vals)
                    m = m + 1
                    call trace_and_pull_back(ccs, rho_vals(ir), theta_vals(it), &
                                             zeta_vals(iz), iota_est(m), drift(m), &
                                             nfail_trace)
                end do
            end do
        end do
        nerrors = nerrors + nfail_trace
        if (nfail_trace == 0) print *, &
            "  PASS: all pullbacks LOCATED on 27 launches x", &
            n_transit*n_per_transit, " samples"

        ! Gate (b): bounded flux-label drift.
        nfail = 0
        if (maxval(drift) > tol_drift) then
            print *, "  FAIL: flux-label drift ", maxval(drift), " > ", tol_drift
            nfail = nfail + 1
        end if
        nerrors = nerrors + nfail
        if (nfail == 0) print *, "  PASS: flux-label drift over ", n_transit, &
            " transits bounded by ", tol_drift, " (max ", maxval(drift), ")"

        ! Gate (c): launch-independent rotational transform.
        nfail = 0
        iota_min = minval(iota_est)
        iota_max = maxval(iota_est)
        if (iota_max - iota_min > tol_iota_spread) then
            print *, "  FAIL: iota ensemble spread ", iota_max - iota_min, &
                " > ", tol_iota_spread
            nfail = nfail + 1
        end if
        if (maxval(abs(iota_est - iota_chart_ref)) > tol_iota_abs) then
            print *, "  FAIL: iota deviates from reference ", iota_chart_ref, &
                " max err ", maxval(abs(iota_est - iota_chart_ref))
            nfail = nfail + 1
        end if
        nerrors = nerrors + nfail
        if (nfail == 0) print *, "  PASS: transform launch-independent, iota=", &
            0.5_dp*(iota_min + iota_max), " ref=", iota_chart_ref, &
            " spread=", iota_max - iota_min
    end subroutine run_ensemble_gates

    ! Trace the closed-form field line launched from chart point (rho0, theta0,
    ! zeta0) over n_transit toroidal transits, pulling every sample back through
    ! the warm chart inverse. Returns the endpoint transform estimate and the
    ! maximum flux-label drift; counts non-LOCATED pullbacks into nfail.
    subroutine trace_and_pull_back(ccs, rho0, theta0, zeta0, iota_out, drift_out, &
                                   nfail)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        real(dp), intent(in) :: rho0, theta0, zeta0
        real(dp), intent(out) :: iota_out, drift_out
        integer, intent(inout) :: nfail

        real(dp) :: x0(3), x(3), u(3), u_guess(3)
        real(dp) :: r, alpha0, phi0, rmaj0
        real(dp) :: alpha, phi, rmaj
        real(dp) :: th_prev, ze_prev, th_un, ze_un, th_start, ze_start
        real(dp) :: dphi, dth, dze
        integer :: k, status, nsteps

        call circular_chart_position(rho0, theta0, zeta0, 1, x0)
        rmaj0 = sqrt(x0(1)**2 + x0(2)**2)
        phi0 = atan2(x0(2), x0(1))
        r = sqrt((rmaj0 - CIRC_R_AXIS)**2 + x0(3)**2)
        alpha0 = atan2(x0(3), rmaj0 - CIRC_R_AXIS)

        nsteps = n_transit*n_per_transit
        dphi = TWOPI/real(n_per_transit, dp)

        u_guess = [rho0 + 0.02_dp, theta0 + 0.05_dp, zeta0 - 0.03_dp]
        drift_out = 0.0_dp
        th_un = 0.0_dp
        ze_un = 0.0_dp
        th_prev = 0.0_dp
        ze_prev = 0.0_dp
        th_start = 0.0_dp
        ze_start = 0.0_dp

        do k = 0, nsteps
            phi = phi0 + real(k, dp)*dphi
            alpha = alpha0 + iota_geom*real(k, dp)*dphi
            rmaj = CIRC_R_AXIS + r*cos(alpha)
            x = [rmaj*cos(phi), rmaj*sin(phi), r*sin(alpha)]

            call ccs%invert_cart(x, u_guess, u, status)
            if (status /= CHARTMAP_LOCATED) then
                print *, "  FAIL: pullback status=", status, " at k=", k, &
                    " launch rho0=", rho0
                nfail = nfail + 1
                iota_out = huge(1.0_dp)
                drift_out = huge(1.0_dp)
                return
            end if
            u_guess = u

            drift_out = max(drift_out, abs(u(1) - rho0))
            if (k == 0) then
                th_un = u(2)
                ze_un = u(3)
                th_start = th_un
                ze_start = ze_un
            else
                dth = modulo(u(2) - th_prev + 0.5_dp*TWOPI, TWOPI) - 0.5_dp*TWOPI
                dze = modulo(u(3) - ze_prev + 0.5_dp*TWOPI, TWOPI) - 0.5_dp*TWOPI
                th_un = th_un + dth
                ze_un = ze_un + dze
            end if
            th_prev = u(2)
            ze_prev = u(3)
        end do

        iota_out = (th_un - th_start)/(ze_un - ze_start)
    end subroutine trace_and_pull_back

end program test_chartmap_gate_fieldline
