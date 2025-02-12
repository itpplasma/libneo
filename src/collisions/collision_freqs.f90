module libneo_collisions

    use libneo_kinds, only : real_kind
    use math_constants, only : pi, ev_to_cgs
    implicit none

    contains

    subroutine calc_coulomb_log()

        implicit none

    end subroutine

    subroutine calc_perp_coll_freq(vel, mass_a, mass_b, charge_a, charge_b, &
        dens_b, temp_b, coulomb_log, coll_freq)

        ! determines the perpendicular collision frequency between two species a and b
        ! a is the test species, b is the background species
        ! units are in cgs (temperature is in eV)

        use stdlib_specialfunctions_gamma, only: lig => lower_incomplete_gamma

        implicit none
        real(kind=real_kind), intent(in) :: vel, mass_a, mass_b, charge_a, charge_b, &
            dens_b, temp_b, coulomb_log
        real(kind=real_kind), intent(out) :: coll_freq

        real(kind=real_kind) :: nue_0, z_ab
        real(kind=real_kind) :: psi_of_x

        real(kind=real_kind) :: p
        p = 1.5d0

        nue_0 = 4.0d0 * pi *charge_a**2 * charge_b**2 * coulomb_log * dens_b / &
            (mass_a**2 * vel**3)
        z_ab = mass_b * vel**2 /(2.0d0 * temp_b * ev_to_cgs)
        
        psi_of_x = lig(p, z_ab)

        coll_freq = nue_0

    end subroutine

    subroutine calc_perp_coll_freq_slow_limit_ee(vel, dens_e, temp_b, coulomb_log, coll_freq)

        use math_constants, only: m_e

        implicit none

        real(kind=real_kind), intent(in) :: vel, dens_e, temp_b, coulomb_log
        real(kind=real_kind), intent(out) :: coll_freq

        coll_freq = 5.8 * 10**(-6d0) * dens_e * coulomb_log * temp_b**-0.5 * (m_e*vel**2 * 0.5 /ev_to_cgs)**(-1.0d0)

    end subroutine

    subroutine calc_perp_coll_freq_fast_limit_ee(vel, dens_e, coulomb_log, coll_freq)

        use math_constants, only: m_e

        implicit none
        real(kind=real_kind), intent(in) :: vel, dens_e, coulomb_log
        real(kind=real_kind), intent(out) :: coll_freq

        coll_freq = 7.7 * 10**(-6d0) * dens_e * coulomb_log * (m_e * vel**2.0d0 * 0.5 / ev_to_cgs)**(-1.5d0)

    end subroutine

end module