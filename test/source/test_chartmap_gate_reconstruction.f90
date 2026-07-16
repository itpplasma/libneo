program test_chartmap_gate_reconstruction
    ! Boozer/SFL chart gate (d) (ROADMAP straight-field-line gauge gates):
    ! converged field reconstruction against the native evaluator under chart
    ! grid refinement.
    !
    ! Pipeline under test, at three chart resolutions (each doubling the grid):
    !   native field sampled on the chart nodes -> periodic batch splines of
    !   (B_R, B_phi, B_Z, |B|) over the chart -> physical probe point x ->
    !   chart inverse invert_cart(x) -> spline evaluation at the inverted
    !   logical point -> compare against the native evaluator at x.
    ! This is the consumer path (GORILLA/SIMPLE style): position inverse plus
    ! field interpolation together, so chart geometry error, inverse residual,
    ! and field spline error all enter the reconstruction error.
    !
    ! Gates: the max relative error (vector B and |B|) decreases monotonically
    ! with refinement, each refinement contracts it by at least 4 (cubic splines
    ! give ~16 asymptotically), and the finest level is below 8e-5 -- the same
    ! bound as the WP3 analytic |B| roundtrip gate.
    !
    ! Scope (honest limits): the native evaluator is the analytic circular
    ! tokamak field on the analytic circular Boozer-convention fixture; the
    ! chart-relevant machinery (node sampling, periodic closure, seams, inverse,
    ! spline order) is production, but there is no JOREK-native comparison here.
    ! libneo has no JOREK-to-chartmap builder yet, so a JOREK-native gate (d)
    ! (against jorek_field on an AUG restart) must follow with that builder.
    ! This test does NOT cover: shaped/X-point geometry, non-nested regions,
    ! fields with toroidal ripple, or the Fourier Boozer-transform pipeline.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinates, only: coordinate_system_t, chartmap_coordinate_system_t, &
                                  make_chartmap_coordinate_system, CHARTMAP_LOCATED
    use math_constants, only: TWOPI
    use interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                           evaluate_batch_splines_3d
    use neo_circular_tokamak_field, only: circular_tokamak_field_t
    use chartmap_gate_fixtures, only: write_circular_boozer_chartmap, &
                                      circular_chart_position, CIRC_R_AXIS
    implicit none

    real(dp), parameter :: B_POL = 0.3183099_dp
    real(dp), parameter :: B_TOR = 1.0_dp
    integer, parameter :: n_levels = 3
    integer, parameter :: nrho_l(n_levels) = [9, 17, 33]
    integer, parameter :: ntheta_l(n_levels) = [16, 32, 64]
    integer, parameter :: nzeta_l(n_levels) = [8, 16, 32]
    real(dp), parameter :: contraction_min = 4.0_dp
    real(dp), parameter :: tol_finest = 8.0e-5_dp

    type(circular_tokamak_field_t) :: field
    real(dp) :: err_vec(n_levels), err_mod(n_levels)
    integer :: level, nerrors

    nerrors = 0
    call field%circular_tokamak_field_init(R_axis=CIRC_R_AXIS, Z_axis=0.0_dp, &
                                           B_pol_ampl=B_POL, B_tor_ampl=B_TOR)

    do level = 1, n_levels
        call reconstruction_error(level, err_vec(level), err_mod(level), nerrors)
        print *, "  level ", level, " grid ", nrho_l(level), ntheta_l(level), &
            nzeta_l(level), " max rel err vec=", err_vec(level), &
            " mod=", err_mod(level)
    end do

    call check_convergence("vector B", err_vec, nerrors)
    call check_convergence("|B|", err_mod, nerrors)

    if (nerrors > 0) then
        print *, "FAILED: ", nerrors, " error(s) detected in reconstruction gate"
        error stop 1
    end if

    print *, "All chartmap field-reconstruction gate tests passed!"

contains

    subroutine check_convergence(label, err, nerrors)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: err(n_levels)
        integer, intent(inout) :: nerrors

        integer :: l, nfail

        nfail = 0
        do l = 1, n_levels - 1
            if (err(l + 1) >= err(l)) then
                print *, "  FAIL: ", label, " error not monotone at level ", l + 1
                nfail = nfail + 1
            else if (err(l)/err(l + 1) < contraction_min) then
                print *, "  FAIL: ", label, " contraction ", err(l)/err(l + 1), &
                    " < ", contraction_min
                nfail = nfail + 1
            end if
        end do
        if (err(n_levels) > tol_finest) then
            print *, "  FAIL: ", label, " finest error ", err(n_levels), &
                " > ", tol_finest
            nfail = nfail + 1
        end if
        nerrors = nerrors + nfail
        if (nfail == 0) print *, "  PASS: ", label, &
            " reconstruction converges under refinement, finest below ", tol_finest
    end subroutine check_convergence

    ! Build the chart and the field splines at one resolution and measure the
    ! max relative reconstruction error over a fixed physical probe set.
    subroutine reconstruction_error(level, emax_vec, emax_mod, nerrors)
        integer, intent(in) :: level
        real(dp), intent(out) :: emax_vec, emax_mod
        integer, intent(inout) :: nerrors

        character(len=64) :: path
        class(coordinate_system_t), allocatable :: cs
        type(BatchSplineData3D) :: spl_b
        real(dp), parameter :: rho_probe(4) = [0.13_dp, 0.38_dp, 0.62_dp, 0.87_dp]
        real(dp) :: u_probe(3), u_guess(3), u_inv(3), x_probe(3)
        real(dp) :: b_true(3), bmod_true, b_recon(4)
        real(dp) :: theta_p, zeta_p, ev, em
        integer :: ip, jp, kp, status, nfail

        write (path, '(a,i0,a)') "chartmap_gate_recon_L", level, ".nc"
        call write_circular_boozer_chartmap(trim(path), nrho_l(level), &
                                            ntheta_l(level), nzeta_l(level), 1)
        call make_chartmap_coordinate_system(cs, trim(path))

        emax_vec = 0.0_dp
        emax_mod = 0.0_dp
        nfail = 0

        select type (ccs => cs)
        type is (chartmap_coordinate_system_t)
            call build_field_splines(nrho_l(level), ntheta_l(level), &
                                     nzeta_l(level), spl_b)

            do kp = 0, 5
                zeta_p = 0.11_dp + 0.97_dp*real(kp, dp)
                do jp = 0, 7
                    theta_p = 0.21_dp + 0.77_dp*real(jp, dp)
                    do ip = 1, size(rho_probe)
                        u_probe = [rho_probe(ip), theta_p, zeta_p]
                        call circular_chart_position(u_probe(1), u_probe(2), &
                                                     u_probe(3), 1, x_probe)
                        ! Inward radial offset keeps the guess out of the warm
                        ! edge band so all probes take the interior 3-D solve.
                        u_guess = [u_probe(1) - 0.05_dp, u_probe(2) + 0.11_dp, &
                                   u_probe(3) - 0.07_dp]
                        call ccs%invert_cart(x_probe, u_guess, u_inv, status)
                        if (status /= CHARTMAP_LOCATED) then
                            print *, "  FAIL: probe not LOCATED, u=", u_probe, &
                                " status=", status
                            nfail = nfail + 1
                            cycle
                        end if
                        call evaluate_batch_splines_3d(spl_b, u_inv, b_recon)
                        call native_field_cyl(x_probe, b_true, bmod_true)
                        ev = sqrt(sum((b_recon(1:3) - b_true)**2))/bmod_true
                        em = abs(b_recon(4) - bmod_true)/bmod_true
                        emax_vec = max(emax_vec, ev)
                        emax_mod = max(emax_mod, em)
                    end do
                end do
            end do
        class default
            print *, "  FAIL: fixture did not load as chartmap type"
            nfail = nfail + 1
        end select
        nerrors = nerrors + nfail
    end subroutine reconstruction_error

    ! Periodic batch splines of (B_R, B_phi, B_Z, |B|) on the chart grid: nodes
    ! at the exact analytic positions, theta/zeta closed with one appended wrap
    ! plane, same order and periodicity layout as the chart splines themselves.
    subroutine build_field_splines(nrho, ntheta, nzeta, spl_b)
        integer, intent(in) :: nrho, ntheta, nzeta
        type(BatchSplineData3D), intent(out) :: spl_b

        real(dp), allocatable :: batch(:, :, :, :)
        real(dp) :: rho, theta, zeta, x(3), b_cyl(3), bmod
        real(dp) :: x_min(3), x_max(3)
        integer :: ir, it, iz

        allocate (batch(nrho, ntheta + 1, nzeta + 1, 4))
        do iz = 1, nzeta
            zeta = TWOPI*real(iz - 1, dp)/real(nzeta, dp)
            do it = 1, ntheta
                theta = TWOPI*real(it - 1, dp)/real(ntheta, dp)
                do ir = 1, nrho
                    rho = real(ir - 1, dp)/real(nrho - 1, dp)
                    call circular_chart_position(rho, theta, zeta, 1, x)
                    call native_field_cyl(x, b_cyl, bmod)
                    batch(ir, it, iz, 1:3) = b_cyl
                    batch(ir, it, iz, 4) = bmod
                end do
            end do
        end do
        batch(:, ntheta + 1, 1:nzeta, :) = batch(:, 1, 1:nzeta, :)
        batch(:, :, nzeta + 1, :) = batch(:, :, 1, :)

        x_min = [0.0_dp, 0.0_dp, 0.0_dp]
        x_max = [1.0_dp, TWOPI, TWOPI]
        call construct_batch_splines_3d(x_min, x_max, batch, [3, 3, 3], &
                                        [.false., .true., .true.], spl_b)
    end subroutine build_field_splines

    ! Native evaluator at a Cartesian point: cylindrical components and modulus
    ! of the production circular tokamak field.
    subroutine native_field_cyl(x, b_cyl, bmod)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b_cyl(3)
        real(dp), intent(out) :: bmod

        real(dp) :: rmaj, phi

        rmaj = sqrt(x(1)**2 + x(2)**2)
        phi = atan2(x(2), x(1))
        call field%compute_bfield([rmaj, phi, x(3)], b_cyl)
        bmod = sqrt(sum(b_cyl**2))
    end subroutine native_field_cyl

end program test_chartmap_gate_reconstruction
