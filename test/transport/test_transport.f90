program test_transport

    use libneo_kinds, only: real_kind
    use libneo_transport, only: calc_D_one_over_nu_11e, init_gauss_laguerre_integration

    implicit none

    call test_D_one_over_nu_11e

    contains

    subroutine test_D_one_over_nu_11e
        
        use libneo_species, only: species_t, init_deuterium_plasma
        use math_constants, only: c, e, ev_to_cgs

        implicit none

        type(species_t) :: species_array(2)
        real(kind=real_kind) :: D11
        real(kind=real_kind) :: R0 = 100.0d0
        real(kind=real_kind) :: B0 = 1.0d4 ! 1 tesla

        call init_gauss_laguerre_integration(5.0d0/2.0d0)

        call init_deuterium_plasma(1.0d3, 1.0d3, 1.0d13, species_array)

        species_array(1)%rho_L = species_array(1)%mass * c * &
            sqrt(species_array(1)%temp * ev_to_cgs / species_array(1)%mass)/ e / B0
        print *, species_array(1)%rho_L

        call calc_D_one_over_nu_11e(2, species_array, R0, D11)
        print *, D11

    end subroutine

end program