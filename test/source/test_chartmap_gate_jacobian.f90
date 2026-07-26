program test_chartmap_gate_jacobian
    ! Boozer/SFL chart gate (a): signed periodic Jacobian positivity across chart
    ! seams on the full torus (ROADMAP straight-field-line gauge gates).
    !
    ! Uses the signed accessor cs%jacobian_det (det of the analytic covariant
    ! basis), not sqrt(abs(det g)), so an orientation flip cannot hide behind the
    ! absolute value. Checks on synthetic analytic fixtures:
    !   - right-handed nfp = 3 rotating-frame chart: det > 0 at every probe of a
    !     full-torus scan that straddles every field-period seam, the theta seam,
    !     and the 2pi zeta wrap;
    !   - seam pairs are continuous (relative jump at solver-noise level, not O(1));
    !   - |jacobian_det| agrees with metric_tensor's sqrtg;
    !   - the mirrored (left-handed) chart gives det < 0 everywhere: the accessor
    !     is genuinely signed;
    !   - right-handed circular nfp = 1 chart: det > 0 including both angular seams.
    ! The magnetic axis rho = 0 is excluded: det -> 0 there by construction of any
    ! polar chart; the scan starts at small but finite rho.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
                                  make_chartmap_coordinate_system
    use math_constants, only: TWOPI
    use chartmap_gate_fixtures, only: write_circular_boozer_chartmap, &
                                      write_stellarator_boozer_chartmap
    implicit none

    integer :: nerrors

    nerrors = 0

    call run_stellarator_checks(nerrors)
    call run_mirrored_checks(nerrors)
    call run_circular_checks(nerrors)

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) detected in Jacobian gate tests"
        error stop 1
    end if

    print *, "All chartmap Jacobian gate tests passed!"

contains

    subroutine load_chart(path, cs)
        character(len=*), intent(in) :: path
        class(coordinate_system_t), allocatable, intent(out) :: cs

        call make_chartmap_coordinate_system(cs, path)
    end subroutine load_chart

    subroutine run_stellarator_checks(nerrors)
        integer, intent(inout) :: nerrors

        character(len=*), parameter :: path = "chartmap_gate_jacobian_rh.nc"
        class(coordinate_system_t), allocatable :: cs

        call write_stellarator_boozer_chartmap(path, 9, 16, 12, 1)
        call load_chart(path, cs)
        select type (ccs => cs)
        type is (chartmap_coordinate_system_t)
            call check_positive_full_torus(ccs, "stellarator RH", nerrors)
            call check_seam_continuity(ccs, "stellarator RH", nerrors)
            call check_sqrtg_consistency(ccs, "stellarator RH", nerrors)
        class default
            print *, "  FAIL: fixture did not load as chartmap type"
            nerrors = nerrors + 1
        end select
    end subroutine run_stellarator_checks

    subroutine run_mirrored_checks(nerrors)
        integer, intent(inout) :: nerrors

        character(len=*), parameter :: path = "chartmap_gate_jacobian_lh.nc"
        class(coordinate_system_t), allocatable :: cs
        real(dp) :: u(3), det
        integer :: ir, it, iz, nfail
        real(dp), parameter :: rho_vals(3) = [0.2_dp, 0.6_dp, 1.0_dp]
        real(dp), parameter :: theta_vals(3) = [0.5_dp, 2.7_dp, 5.0_dp]
        real(dp), parameter :: zeta_fracs(4) = [0.1_dp, 0.35_dp, 0.6_dp, 0.85_dp]

        call write_stellarator_boozer_chartmap(path, 9, 16, 12, -1)
        call load_chart(path, cs)
        nfail = 0
        select type (ccs => cs)
        type is (chartmap_coordinate_system_t)
            do iz = 1, size(zeta_fracs)
                do it = 1, size(theta_vals)
                    do ir = 1, size(rho_vals)
                        u = [rho_vals(ir), theta_vals(it), zeta_fracs(iz)*TWOPI]
                        det = ccs%jacobian_det(u)
                        if (det >= 0.0_dp) then
                            print *, "  FAIL: mirrored chart det>=0 at u=", u, &
                                " det=", det
                            nfail = nfail + 1
                        end if
                    end do
                end do
            end do
        class default
            print *, "  FAIL: fixture did not load as chartmap type"
            nfail = nfail + 1
        end select
        nerrors = nerrors + nfail
        if (nfail == 0) print *, "  PASS: mirrored chart det < 0 (accessor is signed)"
    end subroutine run_mirrored_checks

    subroutine run_circular_checks(nerrors)
        integer, intent(inout) :: nerrors

        character(len=*), parameter :: path = "chartmap_gate_jacobian_circ.nc"
        class(coordinate_system_t), allocatable :: cs

        call write_circular_boozer_chartmap(path, 9, 16, 8, 1)
        call load_chart(path, cs)
        select type (ccs => cs)
        type is (chartmap_coordinate_system_t)
            call check_positive_full_torus(ccs, "circular RH", nerrors)
            call check_seam_continuity(ccs, "circular RH", nerrors)
            call check_sqrtg_consistency(ccs, "circular RH", nerrors)
        class default
            print *, "  FAIL: fixture did not load as chartmap type"
            nerrors = nerrors + 1
        end select
    end subroutine run_circular_checks

    ! Full-torus zeta scan including probes an epsilon on both sides of every
    ! field-period seam and of the 2pi wrap, plus theta-seam probes.
    subroutine check_positive_full_torus(ccs, label, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nerrors

        real(dp), parameter :: eps = 1.0e-6_dp
        real(dp), parameter :: rho_vals(6) = &
            [0.05_dp, 0.25_dp, 0.45_dp, 0.65_dp, 0.85_dp, 1.0_dp]
        real(dp), parameter :: theta_vals(7) = &
            [eps, 0.9_dp, 1.9_dp, 3.1_dp, 4.2_dp, 5.4_dp, TWOPI - eps]
        real(dp) :: zeta_vals(24 + 2*(3 + 1))
        real(dp) :: u(3), det, period
        integer :: ir, it, iz, k, nz, nfail

        period = TWOPI/real(ccs%num_field_periods, dp)
        nz = 0
        do k = 0, 23
            nz = nz + 1
            zeta_vals(nz) = TWOPI*real(k, dp)/24.0_dp + 0.013_dp
        end do
        do k = 0, ccs%num_field_periods
            nz = nz + 1
            zeta_vals(nz) = real(k, dp)*period - eps
            nz = nz + 1
            zeta_vals(nz) = real(k, dp)*period + eps
        end do

        nfail = 0
        do iz = 1, nz
            do it = 1, size(theta_vals)
                do ir = 1, size(rho_vals)
                    u = [rho_vals(ir), theta_vals(it), zeta_vals(iz)]
                    det = ccs%jacobian_det(u)
                    if (det /= det) then
                        print *, "  FAIL: ", label, " det is NaN at u=", u
                        nfail = nfail + 1
                    else if (det <= 0.0_dp) then
                        print *, "  FAIL: ", label, " det<=0 at u=", u, " det=", det
                        nfail = nfail + 1
                    end if
                end do
            end do
        end do
        nerrors = nerrors + nfail
        if (nfail == 0) print *, "  PASS: ", label, &
            " signed Jacobian positive on the full torus incl. seams"
    end subroutine check_positive_full_torus

    ! Seam pairs must agree to solver-noise level: a chart whose orientation or
    ! periodic closure breaks at a seam shows an O(1) jump there.
    subroutine check_seam_continuity(ccs, label, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nerrors

        real(dp), parameter :: eps = 1.0e-6_dp
        real(dp), parameter :: tol_pair = 5.0e-6_dp
        real(dp), parameter :: tol_wrap = 1.0e-10_dp
        real(dp) :: u_lo(3), u_hi(3), det_lo, det_hi, rel
        real(dp) :: period
        integer :: k, nfail

        period = TWOPI/real(ccs%num_field_periods, dp)
        nfail = 0

        ! Field-period seams (interior continuity of the periodic spline).
        do k = 0, ccs%num_field_periods
            u_lo = [0.6_dp, 1.1_dp, real(k, dp)*period - eps]
            u_hi = [0.6_dp, 1.1_dp, real(k, dp)*period + eps]
            det_lo = ccs%jacobian_det(u_lo)
            det_hi = ccs%jacobian_det(u_hi)
            rel = abs(det_hi - det_lo)/max(abs(det_lo), abs(det_hi))
            if (rel > tol_pair) then
                print *, "  FAIL: ", label, " zeta-seam jump at k=", k, " rel=", rel
                nfail = nfail + 1
            end if
        end do

        ! Theta seam.
        u_lo = [0.6_dp, TWOPI - eps, 0.37_dp*period]
        u_hi = [0.6_dp, eps, 0.37_dp*period]
        det_lo = ccs%jacobian_det(u_lo)
        det_hi = ccs%jacobian_det(u_hi)
        rel = abs(det_hi - det_lo)/max(abs(det_lo), abs(det_hi))
        if (rel > tol_pair) then
            print *, "  FAIL: ", label, " theta-seam jump rel=", rel
            nfail = nfail + 1
        end if

        ! Exact 2pi periodicity in both angles (same spline cell, rotation only).
        u_lo = [0.6_dp, 1.3_dp, 0.61_dp*period]
        u_hi = [0.6_dp, 1.3_dp + TWOPI, 0.61_dp*period + TWOPI]
        det_lo = ccs%jacobian_det(u_lo)
        det_hi = ccs%jacobian_det(u_hi)
        rel = abs(det_hi - det_lo)/abs(det_lo)
        if (rel > tol_wrap) then
            print *, "  FAIL: ", label, " 2pi wrap not exact, rel=", rel
            nfail = nfail + 1
        end if

        nerrors = nerrors + nfail
        if (nfail == 0) print *, "  PASS: ", label, " Jacobian continuous across seams"
    end subroutine check_seam_continuity

    ! sqrtg from metric_tensor is sqrt(abs(det g)) built from the same basis;
    ! |jacobian_det| must reproduce it to rounding.
    subroutine check_sqrtg_consistency(ccs, label, nerrors)
        type(chartmap_coordinate_system_t), intent(in) :: ccs
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nerrors

        real(dp), parameter :: tol = 1.0e-9_dp
        real(dp), parameter :: rho_vals(3) = [0.2_dp, 0.55_dp, 0.9_dp]
        real(dp), parameter :: theta_vals(3) = [0.7_dp, 2.9_dp, 5.2_dp]
        real(dp), parameter :: zeta_fracs(3) = [0.2_dp, 0.5_dp, 0.8_dp]
        real(dp) :: u(3), det, g(3, 3), ginv(3, 3), sqrtg, rel
        integer :: ir, it, iz, nfail

        nfail = 0
        do iz = 1, size(zeta_fracs)
            do it = 1, size(theta_vals)
                do ir = 1, size(rho_vals)
                    u = [rho_vals(ir), theta_vals(it), zeta_fracs(iz)*TWOPI]
                    det = ccs%jacobian_det(u)
                    call ccs%metric_tensor(u, g, ginv, sqrtg)
                    rel = abs(abs(det) - sqrtg)/sqrtg
                    if (rel > tol) then
                        print *, "  FAIL: ", label, " |det| vs sqrtg rel=", rel, &
                            " at u=", u
                        nfail = nfail + 1
                    end if
                end do
            end do
        end do
        nerrors = nerrors + nfail
        if (nfail == 0) print *, "  PASS: ", label, " |jacobian_det| matches sqrtg"
    end subroutine check_sqrtg_consistency

end program test_chartmap_gate_jacobian
