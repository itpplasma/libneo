program test_analytical_geoflux
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use analytical_geoflux_field, only: init_analytical_geoflux, &
        splint_analytical_geoflux_field
    implicit none

    real(dp) :: R0, epsilon, kappa, delta, A_param, B0
    integer :: Nripple
    real(dp) :: a0, alpha0, delta0, z0
    real(dp) :: s_tor, theta_geo, phi
    real(dp) :: Acov(3), hcov(3), Bmod, sqgBctr(3)
    integer :: i

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

    print *, 'Initializing analytical geoflux field...'
    call init_analytical_geoflux(R0, epsilon, kappa, delta, A_param, B0, &
                                Nripple, a0, alpha0, delta0, z0)

    print *, 'Testing field evaluation at various locations...'

    s_tor = 0.5d0
    theta_geo = 0.0d0
    phi = 0.0d0

    call splint_analytical_geoflux_field(s_tor, theta_geo, phi, Acov, hcov, &
                                        Bmod, sqgBctr)

    print *, 's=', s_tor, ' theta=', theta_geo, ' phi=', phi
    print *, 'Bmod=', Bmod, ' T'
    print *, 'hcov=', hcov
    print *, 'Acov=', Acov
    print *, 'sqrtg=', sqgBctr(1)

    if (Bmod < 1.0d0 .or. Bmod > 10.0d0) then
        error stop 'Bmod out of expected range'
    end if

    s_tor = 0.9d0
    theta_geo = 1.57d0
    phi = 1.0d0

    call splint_analytical_geoflux_field(s_tor, theta_geo, phi, Acov, hcov, &
                                        Bmod, sqgBctr)

    print *, ''
    print *, 's=', s_tor, ' theta=', theta_geo, ' phi=', phi
    print *, 'Bmod=', Bmod, ' T'
    print *, 'hcov=', hcov
    print *, 'Acov=', Acov
    print *, 'sqrtg=', sqgBctr(1)

    if (Bmod < 1.0d0 .or. Bmod > 10.0d0) then
        error stop 'Bmod out of expected range'
    end if

    print *, ''
    print *, 'Testing with ripple (9 coils)...'

    Nripple = 9
    delta0 = 0.10d0

    call init_analytical_geoflux(R0, epsilon, kappa, delta, A_param, B0, &
                                Nripple, a0, alpha0, delta0, z0)

    s_tor = 0.5d0
    theta_geo = 0.0d0

    do i = 0, 8
        phi = real(i, dp) * 2.0d0 * 3.14159265358979d0 / 9.0d0
        call splint_analytical_geoflux_field(s_tor, theta_geo, phi, Acov, &
                                            hcov, Bmod, sqgBctr)
        print '(A,F8.5,A,F10.6,A)', '  phi=', phi, ' rad, Bmod=', Bmod, ' T'
    end do

    print *, ''
    print *, 'Test PASSED: Analytical geoflux field initialized and evaluated'

end program test_analytical_geoflux
