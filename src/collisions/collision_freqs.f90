module libneo_collisions

    use libneo_kinds, only : real_kind
    use math_constants, only : pi, ev_to_cgs
    implicit none

    contains

    subroutine calc_coulomb_log(interaction_type, temp_a, dens_a, mass_a, charge_num_a, &
        temp_b, dens_b,  mass_b, charge_num_b, coulomb_log)
        ! determines the Coulomb logarithm for a given interaction type
        ! units are in cgs (temperature is in eV)
        ! interaction_type: "ee", "ei", "ie", "ii"
        ! species a is electrons for ee, ei and ie; species b is ions

        use math_constants, only : m_e, m_p

        implicit none

        character(len=2), intent(in) :: interaction_type
        real(kind=real_kind), intent(in) :: temp_a, dens_a, temp_b, dens_b
        real(kind=real_kind), intent(in) :: mass_a, mass_b
        integer, intent(in) :: charge_num_a, charge_num_b
        real(kind=real_kind), intent(out) :: coulomb_log

        select case (trim(interaction_type))
            case ("ee")
                coulomb_log = 23.5 - log(sqrt(dens_a) * temp_a**(-5.0/4.0)) - &
                    sqrt(10**(-5.0) + (log(temp_a) - 2.0d0)**2.0/16.0)
            case ("ei", "ie")
                if (temp_b * m_e / mass_b < temp_a .and. temp_a < 10.0 * charge_num_b**2.0) then
                    coulomb_log = 23.0 - log(sqrt(dens_a) * charge_num_b * temp_a**(-1.5))
                elseif (temp_b * m_e / mass_b < 10.0 * charge_num_b**2.0 .and. 10.0 * charge_num_b**2.0 < temp_a) then
                    coulomb_log = 24.0 - log(sqrt(dens_a) * temp_a**(-1.0))
                elseif (temp_a < temp_b * m_e / mass_b) then
                    coulomb_log = 16.0 - log(sqrt(dens_b) * temp_b**(-1.5) * charge_num_b**2.0 * mass_b/m_p)
                end if
            case ("ii")
                coulomb_log = 23.0 - log(charge_num_a * charge_num_b * (mass_a + mass_b) / (mass_a * temp_b + mass_b * temp_a) &
                    * sqrt(dens_a * charge_num_a**2.0 / temp_a + dens_b * charge_num_b**2.0 / temp_b))
            case default
                print *, "Unknown interaction type"
        end select

        if (coulomb_log < 0.0 .or. coulomb_log > 30.0) then
            stop "Coulomb logarithm has unexpected value"
        end if

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