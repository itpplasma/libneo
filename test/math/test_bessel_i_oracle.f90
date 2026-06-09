program test_bessel_i_oracle
    ! Compares neo_bessel_i against GSL (gsl_sf_bessel_In and
    ! gsl_sf_bessel_In_array) on a log grid over |x| in [1e-300, 700].
    ! Values where either side is below 1e-280 count as matching underflow.

    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
    use libneo_kinds, only: dp
    use neo_bessel_i, only: bessel_in, bessel_in_array
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    interface
        function gsl_set_error_handler_off() bind(c) result(old_handler)
            import :: c_ptr
            type(c_ptr) :: old_handler
        end function gsl_set_error_handler_off

        function gsl_sf_bessel_in(n, x) bind(c, name="gsl_sf_bessel_In") &
            result(value)
            import :: c_int, c_double
            integer(c_int), value :: n
            real(c_double), value :: x
            real(c_double) :: value
        end function gsl_sf_bessel_in

        function gsl_sf_bessel_in_array(nmin, nmax, x, result_array) &
            bind(c, name="gsl_sf_bessel_In_array") result(status)
            import :: c_int, c_double
            integer(c_int), value :: nmin, nmax
            real(c_double), value :: x
            real(c_double) :: result_array(*)
            integer(c_int) :: status
        end function gsl_sf_bessel_in_array
    end interface

    integer, parameter :: n_orders = 8
    integer, parameter :: orders(n_orders) = [0, 1, 2, 5, 10, 50, 100, 200]
    integer, parameter :: n_grid = 200
    real(dp), parameter :: x_min = 1.0d-300
    real(dp), parameter :: x_max = 700.0d0
    ! GSL itself deviates from 40-digit mpmath references by up to ~4.3e-13
    ! (e.g. gsl_sf_bessel_In(1, 4.234e-13) where the exact value is x/2 and
    ! this implementation is bit-accurate), so the agreement tolerance is
    ! 5e-13 rather than the 1e-13 used against mpmath in test_bessel_i.
    real(dp), parameter :: rel_tol = 5.0d-13
    real(dp), parameter :: underflow_tol = 1.0d-280

    type(c_ptr) :: old_handler

    old_handler = gsl_set_error_handler_off()

    call test_scalar_vs_gsl
    call test_array_vs_gsl

contains

    function grid_point(i) result(x)
        integer, intent(in) :: i
        real(dp) :: x

        real(dp) :: u

        u = real(i - 1, dp)/real(n_grid - 1, dp)
        x = x_min*(x_max/x_min)**u
    end function grid_point


    subroutine check(n, x, ours, gsl, worst)
        integer, intent(in) :: n
        real(dp), intent(in) :: x, ours, gsl
        real(dp), intent(inout) :: worst

        real(dp) :: err, tol

        if (abs(ours) < underflow_tol .and. abs(gsl) < underflow_tol) return
        err = abs(ours - gsl)/abs(gsl)
        tol = rel_tol
        ! GSL's I_1 small-argument branch carries O(x) relative error
        ! (measured against mpmath: rel err ~ |x| for |x| < 5e-8, peaking
        ! at 4.2e-8); the exact value there is x/2 to double precision.
        if (n == 1 .and. abs(x) < 1.0d-7) tol = max(rel_tol, 2.0d0*abs(x))
        if (tol == rel_tol) worst = max(worst, err)
        if (err > tol) then
            call print_fail
            print *, "n = ", n, ", x = ", x
            print *, "ours ", ours, ", gsl ", gsl, ", rel err ", err
            stop 1
        end if
    end subroutine check


    subroutine test_scalar_vs_gsl
        real(dp) :: x, ours, gsl, worst
        integer :: i, j

        call print_test("test_scalar_vs_gsl")

        worst = 0.0_dp
        do j = 1, n_orders
            do i = 1, n_grid
                x = grid_point(i)
                ours = bessel_in(orders(j), x)
                gsl = gsl_sf_bessel_in(orders(j), x)
                call check(orders(j), x, ours, gsl, worst)
                ours = bessel_in(orders(j), -x)
                gsl = gsl_sf_bessel_in(orders(j), -x)
                call check(orders(j), -x, ours, gsl, worst)
            end do
        end do
        print *, "    max rel err vs GSL: ", worst

        call print_ok
    end subroutine test_scalar_vs_gsl


    subroutine test_array_vs_gsl
        integer, parameter :: nmax = 200

        real(dp) :: x, ours(0:nmax), gsl(0:nmax), worst
        integer :: i, n, status

        call print_test("test_array_vs_gsl")

        worst = 0.0_dp
        do i = 1, n_grid, 5
            x = grid_point(i)
            call bessel_in_array(nmax, x, ours)
            status = gsl_sf_bessel_in_array(0, nmax, x, gsl)
            if (status /= 0) then
                ! gsl_sf_bessel_In_array zeroes the whole output and
                ! reports underflow once any requested order underflows;
                ! fall back to the scalar GSL oracle for such x.
                do n = 0, nmax
                    gsl(n) = gsl_sf_bessel_in(n, x)
                end do
            end if
            do n = 0, nmax
                call check(n, x, ours(n), gsl(n), worst)
            end do
        end do
        print *, "    max rel err vs GSL: ", worst

        call print_ok
    end subroutine test_array_vs_gsl

end program test_bessel_i_oracle
