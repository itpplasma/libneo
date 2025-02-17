module libneo_transport

    use libneo_kinds, only: real_kind

    implicit none

    logical :: debugging = .false.
    real (kind = real_kind), allocatable, dimension ( : ) :: w
    real (kind = real_kind), allocatable, dimension ( : ) :: x
    logical :: write_to_file = .false.
    integer :: gauss_laguerre_order

    type species
        character(len = 16) :: name
        real(kind = real_kind) :: temp
        real(kind = real_kind) :: dens
        real(kind = real_kind) :: mass
        integer :: charge_num
    end type species

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

    subroutine D_one_over_nu_11e(D11)

        use libneo_collisions, only: calc_perp_coll_freq

        implicit none

        real(kind = real_kind) :: coll_freq, coll_freq_tot
        real(kind = real_kind), intent(out) :: D11
        
        integer :: i
        
        D11 = 0.0d0

        do i = 1, gauss_laguerre_order
            vel = sqrt(2 * x(i) * temp_a / mass_a)
            
            D11 = D11 + w(i) * x(i)**2.0d0
        end do


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