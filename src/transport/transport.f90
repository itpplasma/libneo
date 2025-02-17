module libneo_transport

    use libneo_kinds, only: real_kind

    implicit none

    logical :: debugging = .false.
    real (kind = real_kind), allocatable, dimension ( : ) :: w
    real (kind = real_kind), allocatable, dimension ( : ) :: x
    logical :: write_to_file = .false.
    integer :: gauss_laguerre_order = 64
    real(kind = real_kind) :: eps_eff = 0.025d0

    contains

    subroutine init_gauss_laguerre_integration(alpha)

        ! calculate gauss laguerre weights and abscissas
        ! alpha is the exponent of the monomial in the integration
        ! alpha is 5/2 for D11 and 7/2 for D12

        implicit none
        ! calculate Gauss-Laguerre weights and abscissas
        real(kind = real_kind) :: alpha
        real (kind = real_kind) :: a, b, beta
        
        beta = 0.0D+00 !dummy
        a = 0.0d0
        b = 1.0d0

        allocate(w(gauss_laguerre_order), x(gauss_laguerre_order))

        call calculate_gauss_laguerre_rule(gauss_laguerre_order, alpha, beta, &
            a, b, x, w, write_to_file)

    end subroutine

    subroutine calc_D_one_over_nu_11e(num_species, spec_arr, R0, D11)

        use libneo_collisions, only: calc_perp_coll_freq, calc_coulomb_log
        use libneo_species, only: species_t
        use math_constants, only: ev_to_cgs, pi

        implicit none

        integer, intent(in) :: num_species
        type(species_t), intent(in) :: spec_arr(num_species)
        real(kind = real_kind), intent(in) :: R0
        real(kind = real_kind), intent(out) :: D11

        real(kind = real_kind) :: coll_freq, coll_freq_tot
        real(kind = real_kind), dimension(num_species) :: coll_freq_arr
        real(kind = real_kind) :: coulomb_log, vel
        
        integer :: i, ispecies
        
        D11 = 0.0d0

        do i = 1, gauss_laguerre_order

            vel = sqrt(2.0d0 * x(i) * spec_arr(1)%temp * ev_to_cgs / spec_arr(1)%mass)

            do ispecies = 1, num_species
                if (ispecies .eq. 1) then
                    call calc_coulomb_log("ee", spec_arr(1), spec_arr(ispecies), coulomb_log)
                else
                    call calc_coulomb_log("ei", spec_arr(1), spec_arr(ispecies), coulomb_log)
                end if

                call calc_perp_coll_freq(vel, spec_arr(1), spec_arr(ispecies), coulomb_log, coll_freq)

                coll_freq_arr(ispecies) = coll_freq
            end do

            
            coll_freq_tot = sum(coll_freq_arr)
            print *, "i = ", i, " x(i)= ", x(i), " vel =", vel, " coll_freq(1) = ", coll_freq_arr(1), & 
                " coll_freq(2) = ", coll_freq_arr(2), " coll_freq_tot = ", coll_freq_tot
            D11 = D11 + w(i) / coll_freq_tot
        end do

        print *, "D11 integral = ", D11
        D11 = D11 * 8.0d0 * sqrt(8.0d0) / (9.0d0 * pi**1.5d0) &
            * (spec_arr(1)%temp * ev_to_cgs / spec_arr(1)%mass) &
            * spec_arr(1)%rho_L * eps_eff**1.5d0 / R0**2.0d0

    end subroutine


    subroutine calculate_gauss_laguerre_rule(order, alpha, beta, a, b, x, w, write_to_file)

        implicit none

        integer, intent(in) :: order
        real (kind = real_kind), intent(in) :: alpha, beta, a, b
        logical, intent(in) :: write_to_file
        real (kind = real_kind), dimension(order), intent(out) :: x, w

        integer :: kind
        real (kind = real_kind) r(2)
        real (kind = real_kind) r8_huge
        character ( len = 255 ) :: filename

        if (debugging) then
            write ( *, '(a)' ) ' '
            write ( *, '(a,i8)' ) '  ORDER = ', order
            write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
            write ( *, '(a,g14.6)' ) '  A = ', a
            write ( *, '(a,g14.6)' ) '  B = ', b
        end if

        kind = 5
        call cgqf(order, kind, alpha, beta, a, b, x, w)

        if (write_to_file) then
            filename = 'gauss_laguerre'
            write ( *, '(a)' ) '  FILENAME = "' // trim ( filename ) // '".'
            r(1) = a
            r(2) = r8_huge ( ) 
            call rule_write(order, x, w, r, filename)
        end if

    end subroutine

end module