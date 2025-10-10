module analytical_geoflux_field

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use geoflux_coordinates, only: initialize_analytical_geoflux, geoflux_to_cyl, &
        geoflux_get_axis
    use analytical_tokamak_field, only: analytical_circular_eq_t

    implicit none

    type(analytical_circular_eq_t), save :: stored_eq
    real(dp), save :: axis_R = 0.0_dp
    real(dp), save :: axis_Z = 0.0_dp
    real(dp), save :: psi_axis_val = 0.0_dp
    logical, save :: initialized = .false.

contains

    subroutine init_analytical_geoflux(R0, epsilon, kappa, delta, A_param, B0, &
                                      Nripple, a0, alpha0, delta0, z0, &
                                      ns_cache, ntheta_cache)
        real(dp), intent(in) :: R0, epsilon, kappa, delta, A_param, B0
        integer, intent(in) :: Nripple
        real(dp), intent(in) :: a0, alpha0, delta0, z0
        integer, intent(in), optional :: ns_cache, ntheta_cache

        if (initialized) then
            call stored_eq%cleanup()
        end if

        call stored_eq%init(R0, epsilon, kappa, delta, A_param, B0, &
                           Nripple_in=Nripple, a0_in=a0, alpha0_in=alpha0, &
                           delta0_in=delta0, z0_in=z0)

        call initialize_analytical_geoflux(R0, epsilon, kappa, delta, A_param, B0, &
                                          Nripple, a0, alpha0, delta0, z0, &
                                          ns_cache, ntheta_cache)

        call geoflux_get_axis(axis_R, axis_Z)
        psi_axis_val = stored_eq%eval_psi(axis_R, axis_Z)

        initialized = .true.
    end subroutine init_analytical_geoflux

    subroutine splint_analytical_geoflux_field(s_tor, theta_geo, phi, Acov, hcov, &
                                              Bmod, sqgBctr)
        real(dp), intent(in) :: s_tor, theta_geo, phi
        real(dp), intent(out) :: Acov(3), hcov(3), Bmod
        real(dp), intent(out), optional :: sqgBctr(3)

        real(dp) :: x_geo(3)
        real(dp) :: x_cyl(3)
        real(dp) :: jac(3, 3)
        real(dp) :: dRds, dRdt, dZds, dZdt
        real(dp) :: det_s_theta, sqrtg
        real(dp) :: Br, Bphi, Bz
        real(dp) :: Bcov_s, Bcov_theta, Bcov_phi
        real(dp) :: psi_local

        if (.not. initialized) then
            error stop "analytical_geoflux_field: call init_analytical_geoflux first"
        end if

        x_geo = [s_tor, theta_geo, phi]
        call geoflux_to_cyl(x_geo, x_cyl, jac)

        call stored_eq%eval_bfield_ripple(x_cyl(1), phi, x_cyl(3), Br, Bz, Bphi, Bmod)

        psi_local = stored_eq%eval_psi(x_cyl(1), x_cyl(3))

        dRds = jac(1, 1)
        dRdt = jac(1, 2)
        dZds = jac(3, 1)
        dZdt = jac(3, 2)

        det_s_theta = dRds*dZdt - dRdt*dZds
        sqrtg = abs(x_cyl(1)*det_s_theta)

        Bcov_s = Br*dRds + Bz*dZds
        Bcov_theta = Br*dRdt + Bz*dZdt
        Bcov_phi = Bphi*x_cyl(1)

        if (Bmod > 0.0_dp) then
            hcov(1) = Bcov_s/Bmod
            hcov(2) = Bcov_theta/Bmod
            hcov(3) = Bcov_phi/Bmod
        else
            hcov = 0.0_dp
        end if

        Acov(1) = 0.0_dp
        Acov(2) = 0.0_dp
        Acov(3) = psi_local - psi_axis_val

        if (present(sqgBctr)) then
            sqgBctr = 0.0_dp
            sqgBctr(1) = sqrtg
        end if
    end subroutine splint_analytical_geoflux_field

end module analytical_geoflux_field
