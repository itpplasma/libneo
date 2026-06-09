program test_incomplete_gamma_gsl_oracle

    ! compares the in-tree lower_incomplete_gamma against GSL; built only
    ! when GSL is available, so GSL stays a test oracle, not a dependency

    use, intrinsic :: iso_c_binding, only: c_double
    use libneo_kinds, only: dp
    use libneo_collisions, only: lower_incomplete_gamma
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    interface
        function gsl_sf_gamma(x) bind(c, name='gsl_sf_gamma')
            import :: c_double
            implicit none
            real(c_double), value :: x
            real(c_double) :: gsl_sf_gamma
        end function gsl_sf_gamma
    end interface

    interface
        function gsl_sf_gamma_inc_P(a, x) bind(c, name='gsl_sf_gamma_inc_P')
            import :: c_double
            implicit none
            real(c_double), value :: a, x
            real(c_double) :: gsl_sf_gamma_inc_P
        end function gsl_sf_gamma_inc_P
    end interface

    integer, parameter :: n_a = 24, n_x = 81
    real(dp), parameter :: rel_tol = 1.0d-12
    real(dp) :: a, x, ours, oracle
    integer :: i, j

    call print_test("test_incomplete_gamma_gsl_oracle")

    do i = 1, n_a
        a = 0.25d0*i
        do j = 1, n_x
            x = 10.0d0**(-8.0d0 + 0.125d0*(j - 1))
            ours = lower_incomplete_gamma(a, x)
            oracle = gsl_sf_gamma_inc_P(a, x)*gsl_sf_gamma(a)
            if (abs(ours - oracle) > rel_tol*oracle) then
                call print_fail
                print *, "a = ", a, ", x = ", x
                print *, "got ", ours, ", GSL gives ", oracle
                stop "lower_incomplete_gamma deviates from GSL"
            end if
        end do
    end do
    call print_ok

end program test_incomplete_gamma_gsl_oracle
