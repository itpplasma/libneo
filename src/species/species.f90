module libneo_species

    use libneo_kinds, only: real_kind

    implicit none

    type species_t
        character(len=16) :: name ! whatever, e.g. deuterium
        character(len=1) :: typ ! either e or i
        real(kind=real_kind) :: temp
        real(kind=real_kind) :: dens
        real(kind=real_kind) :: mass
        integer :: charge_num
        real(kind=real_kind) :: rho_L
        real(kind=real_kind), allocatable, dimension(:) :: coulomb_log
    end type species_t

    contains

    subroutine init_electron_species(temp, dens, electron)

        use math_constants, only: m_e
        implicit none
        real(kind=real_kind), intent(in) :: temp, dens
        type(species_t), intent(out) :: electron

        electron%name = "electron"
        electron%typ = "e"
        electron%temp = temp
        electron%dens = dens
        electron%mass = m_e
        electron%charge_num = -1

    end subroutine

    subroutine init_deuterium_species(temp, dens, deuterium)

        use math_constants, only: m_D
        implicit none
        real(kind=real_kind), intent(in) :: temp, dens
        type(species_t), intent(out) :: deuterium

        deuterium%name = "deuterium"
        deuterium%typ = "i"
        deuterium%temp = temp
        deuterium%dens = dens
        deuterium%mass = m_D
        deuterium%charge_num = 1

    end subroutine

    subroutine init_deuterium_plasma(temp_e, temp_D, dens, deuterium_plasma)

        implicit none
        real(kind=real_kind), intent(in) :: temp_e, temp_D, dens
        type(species_t), intent(out) :: deuterium_plasma(2)

        call init_electron_species(temp_e, dens, deuterium_plasma(1))
        call init_deuterium_species(temp_D, dens, deuterium_plasma(2))

    end subroutine

end module