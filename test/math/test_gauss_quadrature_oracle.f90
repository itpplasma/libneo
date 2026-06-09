program test_gauss_quadrature_oracle

    ! Compares gauss_legendre_ab against GSL gsl_integration_glfixed_table.
    ! Link with -lgsl (and -lgslcblas where required).

    use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_size_t, c_int, &
                                           c_associated
    use libneo_kinds, only: dp
    use neo_gauss_quadrature, only: gauss_legendre_ab
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    interface
        function gsl_glfixed_table_alloc(n) result(t) &
            bind(c, name="gsl_integration_glfixed_table_alloc")
            import :: c_ptr, c_size_t
            integer(c_size_t), value :: n
            type(c_ptr) :: t
        end function gsl_glfixed_table_alloc

        function gsl_glfixed_point(a, b, i, xi, wi, t) result(status) &
            bind(c, name="gsl_integration_glfixed_point")
            import :: c_ptr, c_size_t, c_double, c_int
            real(c_double), value :: a, b
            integer(c_size_t), value :: i
            real(c_double), intent(out) :: xi, wi
            type(c_ptr), value :: t
            integer(c_int) :: status
        end function gsl_glfixed_point

        subroutine gsl_glfixed_table_free(t) &
            bind(c, name="gsl_integration_glfixed_table_free")
            import :: c_ptr
            type(c_ptr), value :: t
        end subroutine gsl_glfixed_table_free
    end interface

    ! n = 33 and 50 exercise GSL's Newton-iterated (non-precomputed) tables.
    ! Those carry up to ~1e-10 relative weight error (checked against mpmath
    ! at 40 digits, where our weights agree to ~1e-16), hence the looser
    ! tolerance for them
    integer, parameter :: n_list(7) = [4, 8, 16, 32, 33, 50, 64]
    real(dp), parameter :: tol_list(7) = [1.0d-13, 1.0d-13, 1.0d-13, &
                                          1.0d-13, 1.0d-9, 1.0d-9, 1.0d-13]
    real(dp), parameter :: a_list(2) = [-1.0_dp, 0.25_dp]
    real(dp), parameter :: b_list(2) = [1.0_dp, 3.75_dp]

    integer :: j, k
    real(dp) :: max_err_strict, max_err_loose

    max_err_strict = 0.0_dp
    max_err_loose = 0.0_dp
    call print_test("test_gauss_legendre_vs_gsl_glfixed")
    do j = 1, size(n_list)
        do k = 1, size(a_list)
            if (tol_list(j) <= 1.0d-13) then
                call compare_with_gsl(n_list(j), a_list(k), b_list(k), &
                                      tol_list(j), max_err_strict)
            else
                call compare_with_gsl(n_list(j), a_list(k), b_list(k), &
                                      tol_list(j), max_err_loose)
            end if
        end do
    end do
    print *, "    max rel err vs GSL precomputed tables:", max_err_strict
    print *, "    max rel err vs GSL iterated tables:   ", max_err_loose
    call print_ok

contains

    subroutine compare_with_gsl(n, a, b, rel_tol, max_rel_err)
        integer, intent(in) :: n
        real(dp), intent(in) :: a, b, rel_tol
        real(dp), intent(inout) :: max_rel_err

        type(c_ptr) :: table
        real(dp) :: x(n), w(n), xi, wi, err_x, err_w
        integer(c_int) :: status
        integer :: i

        table = gsl_glfixed_table_alloc(int(n, c_size_t))
        if (.not. c_associated(table)) then
            call print_fail
            print *, "GSL table alloc failed for n=", n
            stop 1
        end if
        call gauss_legendre_ab(n, a, b, x, w)
        do i = 1, n
            status = gsl_glfixed_point(a, b, int(i - 1, c_size_t), xi, wi, &
                                       table)
            if (status /= 0) then
                call print_fail
                print *, "GSL glfixed_point failed, n=", n, " i=", i
                stop 1
            end if
            ! nodes can pass through zero on [-1,1]: scale by interval width
            err_x = abs(x(i) - xi)/(b - a)
            err_w = abs(w(i)/wi - 1.0_dp)
            max_rel_err = max(max_rel_err, err_x, err_w)
            if (err_x > rel_tol .or. err_w > rel_tol) then
                call print_fail
                print *, "n=", n, " i=", i, " ours:", x(i), w(i), &
                    " gsl:", xi, wi
                stop 1
            end if
        end do
        call gsl_glfixed_table_free(table)
    end subroutine compare_with_gsl

end program test_gauss_quadrature_oracle
