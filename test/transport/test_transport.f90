program test_transport

    use libneo_kinds, only: real_kind
    use libneo_transport, only: calc_D_one_over_nu_11e, init_gauss_laguerre_integration
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    call test_D_one_over_nu_11e

    contains

    subroutine test_D_one_over_nu_11e
        
        use libneo_species, only: species_t, init_deuterium_plasma
        use math_constants, only: c, e, ev_to_cgs

        implicit none

        type(species_t) :: species_array(2)
        real(kind=real_kind) :: D11
        real(kind=real_kind) :: R0 = 550.0d0 ! values for W-7X
        real(kind=real_kind) :: B0 = 2.5d4 

        call init_gauss_laguerre_integration(5.0d0/2.0d0)
        call init_deuterium_plasma(3.0d3, 1.5d3, 1.0d14, species_array)

        species_array(1)%rho_L = species_array(1)%mass * c * &
            sqrt(species_array(1)%temp * ev_to_cgs / species_array(1)%mass)/ e / B0

        call calc_D_one_over_nu_11e(2, species_array, R0, D11)
        print *, "D11e = ", D11, " compared to ", 0.00049601
        if (abs(D11 - 0.00049601) > 1e-6) then ! number is for the parameters above
            call print_fail
            stop "D11e != 0.00049601"
        else
            call print_ok
        end if

    end subroutine

end program