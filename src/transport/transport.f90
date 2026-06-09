module libneo_transport

    implicit none

    integer, parameter :: dp = kind(1.0d0)

    integer :: gauss_laguerre_order = 32
    real(dp) :: eps_eff = 0.025d0 ! value for W-7X high mirror configuration

    contains

    subroutine init_gauss_laguerre_integration(alpha, weights, x_abscissas)

        use neo_gauss_quadrature, only: gauss_gen_laguerre

        ! calculate gauss laguerre weights and abscissas
        ! alpha is the exponent of the monomial in the integration
        ! alpha is 5/2 for D11 and 7/2 for D12
        ! weights and x_abscissas are the arrays to be filled with the weights and abscissas

        implicit none
        real(dp), intent(in) :: alpha
        real(dp), dimension(gauss_laguerre_order), intent(out) :: weights, x_abscissas

        call gauss_gen_laguerre(gauss_laguerre_order, alpha, x_abscissas, weights)

    end subroutine

    subroutine calc_D_one_over_nu(ind_species, num_species, spec_arr, R0, w, x, D_one_over_nu)

        use libneo_collisions, only: calc_perp_coll_freq
        use libneo_species, only: species_t
        use math_constants, only: ev_to_cgs, pi

        implicit none

        integer, intent(in) :: num_species, ind_species
        type(species_t), intent(in) :: spec_arr(num_species)
        real(dp), intent(in) :: R0
        real(dp), dimension(gauss_laguerre_order), intent(in) :: w, x
        real(dp), intent(out) :: D_one_over_nu

        real(dp) :: coll_freq, coll_freq_tot
        real(dp), dimension(num_species) :: coll_freq_arr
        real(dp) :: vel
        
        integer :: i, ispecies
        
        D_one_over_nu = 0.0d0

        do i = 1, gauss_laguerre_order
            do ispecies = 1, num_species
                vel = sqrt(2.0d0 * x(i) * spec_arr(ispecies)%temp * ev_to_cgs / spec_arr(ispecies)%mass)
                call calc_perp_coll_freq(vel, spec_arr(ind_species), spec_arr(ispecies), &
                    spec_arr(ind_species)%coulomb_log(ispecies), coll_freq)
                coll_freq_arr(ispecies) = coll_freq
            end do
            coll_freq_tot = sum(coll_freq_arr)
            D_one_over_nu = D_one_over_nu + w(i) / coll_freq_tot
        end do

        D_one_over_nu = D_one_over_nu * 8.0d0 * sqrt(8.0d0) / (9.0d0 * pi**1.5d0) &
            * (spec_arr(ind_species)%temp * ev_to_cgs / spec_arr(ind_species)%mass) &
            * spec_arr(ind_species)%rho_L**2.0d0 * eps_eff**1.5d0 / R0**2.0d0

    end subroutine

end module
