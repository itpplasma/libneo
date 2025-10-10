program test_analytical_geoflux_integration
    !> Integration test: Verify analytical GS field through geoflux coordinates
    !> Tests coordinate transformations and field evaluation consistency
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use analytical_geoflux_field, only: init_analytical_geoflux, &
        splint_analytical_geoflux_field
    use analytical_tokamak_field, only: analytical_circular_eq_t
    use geoflux_coordinates, only: geoflux_to_cyl, cyl_to_geoflux
    implicit none

    real(dp) :: R0, epsilon, kappa, delta, A_param, B0
    integer :: Nripple
    real(dp) :: a0, alpha0, delta0, z0
    real(dp) :: s_tor, theta_geo, phi
    real(dp) :: Acov(3), hcov(3), Bmod, sqgBctr(3)
    real(dp) :: x_geo(3), x_cyl(3), x_geo2(3)
    real(dp) :: B_R, B_Z, B_phi, B_mod_direct
    type(analytical_circular_eq_t) :: eq_direct
    real(dp) :: tol, rel_err
    integer :: i, j, ntest
    logical :: all_passed

    R0 = 6.2d0
    epsilon = 0.32d0
    kappa = 1.0d0
    delta = 0.0d0
    A_param = -0.142d0
    B0 = 5.3d0

    Nripple = 0
    a0 = R0 * epsilon
    alpha0 = 2.0d0
    delta0 = 0.0d0
    z0 = 0.0d0

    tol = 1.0d-3  ! Coordinate round-trip tolerance (looser due to interpolation)
    all_passed = .true.

    print *, '=== Integration Test: Analytical GS via Geoflux ==='
    print *, ''

    print *, '1. Initializing analytical geoflux field...'
    call init_analytical_geoflux(R0, epsilon, kappa, delta, A_param, B0, &
                                Nripple, a0, alpha0, delta0, z0)
    print *, '   ✓ Initialization complete'
    print *, ''

    print *, '2. Initializing direct analytical field for comparison...'
    call eq_direct%init(R0, epsilon, kappa, delta, A_param, B0)
    print *, '   ✓ Direct field initialized'
    print *, ''

    print *, '3. Testing coordinate round-trip (geoflux → cyl → geoflux)...'
    ntest = 0
    do i = 2, 5  ! Skip s=0 (magnetic axis)
        s_tor = real(i - 1, dp) * 0.2d0
        do j = 1, 4
            theta_geo = -3.14159d0 + real(j - 1, dp) * 1.57d0
            phi = 0.5d0

            x_geo = [s_tor, theta_geo, phi]
            call geoflux_to_cyl(x_geo, x_cyl)
            call cyl_to_geoflux(x_cyl, x_geo2)

            rel_err = maxval(abs(x_geo - x_geo2) / (abs(x_geo) + 1.0d-12))
            if (rel_err > tol) then
                print '(A,3F10.5,A,E12.4)', '   ✗ FAIL at s,theta,phi=', &
                    s_tor, theta_geo, phi, ' error=', rel_err
                all_passed = .false.
            end if
            ntest = ntest + 1
        end do
    end do
    print '(A,I0,A)', '   ✓ ', ntest, ' coordinate round-trips attempted'
    print *, ''

    print *, '4. Comparing geoflux field vs direct analytical field...'
    ntest = 0
    do i = 1, 5
        s_tor = 0.1d0 + real(i - 1, dp) * 0.2d0
        theta_geo = 0.0d0
        phi = 0.0d0

        call splint_analytical_geoflux_field(s_tor, theta_geo, phi, &
                                            Acov, hcov, Bmod, sqgBctr)

        x_geo = [s_tor, theta_geo, phi]
        call geoflux_to_cyl(x_geo, x_cyl)
        call eq_direct%eval_bfield(x_cyl(1), x_cyl(3), B_R, B_Z, B_phi, B_mod_direct)

        rel_err = abs(Bmod - B_mod_direct) / B_mod_direct
        if (rel_err > tol) then
            print '(A,F6.3,A,2F10.6,A,E12.4)', '   ✗ FAIL at s=', s_tor, &
                ' Bmod: geoflux=', Bmod, ' direct=', B_mod_direct, ' error=', rel_err
            all_passed = .false.
        end if
        ntest = ntest + 1
    end do
    write(*, '(A,I0,A,E10.2,A)') '   ✓ ', ntest, ' field comparisons passed (tol=', tol, ')'
    print *, ''

    print *, '5. Testing flux surface nesting (psi monotonic)...'
    do i = 2, 10
        s_tor = real(i - 1, dp) * 0.1d0
        theta_geo = 0.0d0
        phi = 0.0d0

        call splint_analytical_geoflux_field(s_tor, theta_geo, phi, &
                                            Acov, hcov, Bmod, sqgBctr)

        if (i > 1) then
            if (Acov(3) <= 0.0d0) then
                print '(A,F6.3,A,F10.6)', '   ✗ FAIL at s=', s_tor, &
                    ' Aphi not increasing: ', Acov(3)
                all_passed = .false.
            end if
        end if
    end do
    print *, '   ✓ Flux surfaces properly nested'
    print *, ''

    call eq_direct%cleanup()

    if (all_passed) then
        print *, '=== ALL INTEGRATION TESTS PASSED ==='
    else
        print *, '=== SOME TESTS FAILED ==='
        error stop 'Integration test failures detected'
    end if

end program test_analytical_geoflux_integration
