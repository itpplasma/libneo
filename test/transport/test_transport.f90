program test_transport

    use libneo_kinds, only: dp
    use libneo_transport, only: init_gauss_laguerre_integration, calc_D_one_over_nu, gauss_laguerre_order
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    call test_D_one_over_nu_11e

    contains

    subroutine test_D_one_over_nu_11e
        
        use libneo_species, only: species_t, init_deuterium_plasma
        use libneo_collisions, only: fill_species_arr_coulomb_log
        use math_constants, only: c, e, ev_to_cgs

        implicit none

        integer, parameter :: num_species = 2
        integer :: i
        type(species_t) :: species_array(num_species)
        real(dp) :: D11
        real(dp) :: R0 = 550.0d0 ! values for W-7X
        real(dp) :: B0 = 2.5d4 
        real(dp), allocatable, dimension(:) :: w, x

        call init_deuterium_plasma(3.0d3, 1.5d3, 1.0d14, species_array)

        do i = 1, num_species
            species_array(i)%rho_L = species_array(i)%mass * c * &
                sqrt(species_array(i)%temp * ev_to_cgs / species_array(i)%mass)/ e / B0
        end do
        
        call fill_species_arr_coulomb_log(2, species_array)

        if (.not. allocated(w)) allocate(w(gauss_laguerre_order), x(gauss_laguerre_order))

        call init_gauss_laguerre_integration(5.0d0/2.0d0, w, x)
        call calc_D_one_over_nu(1, 2, species_array, R0, w, x, D11)

        print *, "D11e = ", D11, " compared to ", 0.00049601
        if (abs(1.0d0 - 0.00049601/D11) > 0.05d0) then ! number is for the parameters above
            call print_fail
            stop "D11e != 0.00049601, relative error larger than 5%"
        else
            call print_ok
        end if
        
        call calc_D_one_over_nu(2, 2, species_array, R0, w, x, D11)
        print *, "D11i = ", D11, " compared to ", 1512.81
        if (abs(1.0d0 - 1512.81/D11) > 0.05d0) then ! number is for the parameters above
            call print_fail
            stop "D11e != 1512.81, relative error larger than 5%"
        else
            call print_ok
        end if

        call init_gauss_laguerre_integration(7.0d0/2.0d0, w, x)
        call calc_D_one_over_nu(1, 2, species_array, R0, w, x, D11)
        print *, "D12e = ", D11, " compared to ", 0.00241398
        if (abs(1 - 0.00241398/D11) > 0.05d0) then ! number is for the parameters above
            call print_fail
            stop "D12e != 0.00241398, relative error larger than 5%"
        else
            call print_ok
        end if
        
        call calc_D_one_over_nu(2, 2, species_array, R0, w, x, D11)
        print *, "D12i = ", D11, " compared to ", 7362.44
        if (abs(1.0d0 - 7362.44/D11) > 0.05d0) then ! number is for the parameters above
            call print_fail
            stop "D11e != 7362.44, relative error larger than 5%"
        else
            call print_ok
        end if

    end subroutine

end program