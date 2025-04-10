module libneo_species

    implicit none

    integer, parameter :: dp = kind(1.0d0)

    type species_t
        character(len=16) :: name ! whatever, e.g. deuterium
        character(len=1) :: typ ! either e or i
        integer :: typ_int ! type integer, either =1 for electrons or =2 for ions
        real(dp) :: temp
        real(dp) :: dens
        real(dp) :: mass
        integer :: charge_num
        real(dp) :: rho_L
        real(dp), allocatable, dimension(:) :: coulomb_log
    end type species_t

    contains

    subroutine init_electron_species(temp, dens, electron)

        use math_constants, only: m_e
        implicit none
        real(dp), intent(in) :: temp, dens
        type(species_t), intent(out) :: electron

        electron%name = "electron"
        electron%typ = "e"
        electron%typ_int = 1
        electron%temp = temp
        electron%dens = dens
        electron%mass = m_e
        electron%charge_num = -1

    end subroutine

    subroutine init_deuterium_species(temp, dens, deuterium)

        use math_constants, only: m_D
        implicit none
        real(dp), intent(in) :: temp, dens
        type(species_t), intent(out) :: deuterium

        deuterium%name = "deuterium"
        deuterium%typ = "i"
        deuterium%typ_int = 2
        deuterium%temp = temp
        deuterium%dens = dens
        deuterium%mass = m_D
        deuterium%charge_num = 1

    end subroutine

    subroutine init_deuterium_plasma(temp_e, temp_D, dens, deuterium_plasma)

        implicit none
        real(dp), intent(in) :: temp_e, temp_D, dens
        type(species_t), intent(out) :: deuterium_plasma(2)

        call init_electron_species(temp_e, dens, deuterium_plasma(1))
        call init_deuterium_species(temp_D, dens, deuterium_plasma(2))

    end subroutine

end module
