module libneo_collisions

    use libneo_kinds, only : real_kind
    implicit none

    contains

    subroutine calc_coulomb_log(interaction_type, species_a, species_b, coulomb_log)
        ! determines the Coulomb logarithm for a given interaction type
        ! units are in cgs (temperature is in eV)
        ! interaction_type: "ee", "ei", "ie", "ii"
        ! species a is electrons for ee, ei and ie; species b is ions

        use math_constants, only : m_e, m_p
        use libneo_species, only: species_t

        implicit none

        character(len=2), intent(in) :: interaction_type
        type(species_t), intent(in) :: species_a, species_b
        real(kind=real_kind), intent(out) :: coulomb_log

        select case (trim(interaction_type))
            case ("ee")
                coulomb_log = 23.5 - log(sqrt(species_a%dens) * species_a%temp**(-5.0/4.0)) - &
                    sqrt(10**(-5.0) + (log(species_a%temp) - 2.0d0)**2.0/16.0)
            case ("ei", "ie")
                if (species_b%temp * m_e / species_b%mass < species_a%temp &
                    .and. species_a%temp < 10.0 * species_b%charge_num**2.0) then
                    coulomb_log = 23.0 - log(sqrt(species_a%dens) * species_b%charge_num * &
                    species_a%temp**(-1.5))
                elseif (species_b%temp * m_e / species_b%mass < 10.0 * species_b%charge_num**2.0 &
                    .and. 10.0 * species_b%charge_num**2.0 < species_a%temp) then
                    coulomb_log = 24.0 - log(sqrt(species_a%dens) * species_a%temp**(-1.0))
                elseif (species_a%temp < species_b%temp * m_e / species_b%mass) then
                    coulomb_log = 16.0 - log(sqrt(species_b%dens) * species_b%temp**(-1.5) &
                        * species_b%charge_num**2.0 * species_b%mass/m_p)
                end if
            case ("ii")
                coulomb_log = 23.0 - log(species_a%charge_num * species_b%charge_num &
                    * (species_a%mass + species_b%mass) / (species_a%mass * species_b%temp &
                    + species_b%mass * species_a%temp) &
                    * sqrt(species_a%dens * species_a%charge_num**2.0 / species_a%temp &
                    + species_b%dens * species_b%charge_num**2.0 / species_b%temp))
            case default
                print *, "Unknown interaction type"
        end select

        if (coulomb_log < 0.0 .or. coulomb_log > 30.0) then
            stop "Coulomb logarithm has unexpected value"
        end if

    end subroutine

    subroutine calc_perp_coll_freq(vel, species_a, species_b, coulomb_log, coll_freq)

        ! determines the perpendicular collision frequency between two species a and b
        ! a is the test species, b is the background species
        ! units are in cgs (temperature is in eV)

        use stdlib_specialfunctions_gamma, only: lig => lower_incomplete_gamma
        use math_constants, only: pi, ev_to_cgs, E
        use libneo_species, only: species_t

        implicit none

        real(kind=real_kind), intent(in) :: vel, coulomb_log
        type(species_t), intent(in) :: species_a, species_b
        real(kind=real_kind), intent(out) :: coll_freq
        real(kind=real_kind) :: nue_0, z_ab
        real(kind=real_kind) :: psi_of_x
        real(kind=real_kind) :: p = 1.5d0

        nue_0 = 4.0d0 * pi *species_a%charge_num**2d0 * species_b%charge_num**2d0 * E**4d0 &
            * coulomb_log * species_b%dens / (species_a%mass**2d0 * vel**3d0)
        z_ab = species_b%mass * vel**2d0 /(2.0d0 * species_b%temp * ev_to_cgs)
        psi_of_x = lig(p, z_ab)

        coll_freq = 2d0 * nue_0 * ((1.0d0 - 1.0d0/(2.0d0 * z_ab)) * psi_of_x *2.0d0 /sqrt(pi) &
            + 2.0d0 / sqrt(pi) * z_ab**0.5d0 * exp(-z_ab))

    end subroutine

    subroutine calc_perp_coll_freq_slow_limit_ee(vel, dens_e, temp_b, coulomb_log, coll_freq)

        use math_constants, only: m_e, ev_to_cgs

        implicit none

        real(kind=real_kind), intent(in) :: vel, dens_e, temp_b, coulomb_log
        real(kind=real_kind), intent(out) :: coll_freq

        coll_freq = 5.8d0 * 10d0**(-6d0) * dens_e * coulomb_log * temp_b**-0.5d0 * (m_e*vel**2d0 * 0.5d0 /ev_to_cgs)**(-1.0d0)

    end subroutine

    subroutine calc_perp_coll_freq_fast_limit_ee(vel, dens_e, coulomb_log, coll_freq)

        use math_constants, only: m_e, ev_to_cgs

        implicit none
        real(kind=real_kind), intent(in) :: vel, dens_e, coulomb_log
        real(kind=real_kind), intent(out) :: coll_freq

        coll_freq = 7.7d0 * 10.0d0**(-6.0d0) * dens_e * coulomb_log * (m_e * vel**2.0d0 * 0.5d0 / ev_to_cgs)**(-1.5d0)

    end subroutine

end module