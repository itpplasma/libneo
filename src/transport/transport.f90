module libneo_transport

    use libneo_kinds, only: real_kind

    implicit none

    logical :: debugging = .false.

    contains

    subroutine D_one_over_nu_11()

        use libneo_collisions, only: calc_coulomb_log, calc_perp_coll_freq

        implicit none

        ! calculate Gauss-Laguerre weights and abscissas
        real (kind = real_kind) :: a, alpha, b, beta
        integer :: order
        real ( kind = real_kind ), allocatable, dimension ( : ) :: w
        real ( kind = real_kind ), allocatable, dimension ( : ) :: x
        logical :: write_to_file
        
        beta = 0.0D+00 !dummy
        order = 100
        alpha = 5.0d0/2.0d0
        a = 0.0d0
        b = 1.0d0
        write_to_file = .false.

        allocate(w(order), x(order))

        call calculate_gauss_laguerre_rule(order, alpha, beta, a, b, x, w, write_to_file)

        deallocate(w, x)

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