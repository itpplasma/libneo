module spline1d_many_omp_target
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: spline1d_many_omp_target_setup
    public :: spline1d_many_omp_target_teardown
    public :: spline1d_many_omp_target_eval_resident

    logical :: is_setup = .false.

contains

    subroutine spline1d_many_omp_target_setup(coeff, x, y)
        real(dp), intent(in) :: coeff(:, :, :)
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(:)

        if (is_setup) return
!$omp target enter data map(to: coeff)
!$omp target enter data map(to: x)
!$omp target enter data map(alloc: y)
        is_setup = .true.
    end subroutine spline1d_many_omp_target_setup

    subroutine spline1d_many_omp_target_teardown(coeff, x, y)
        real(dp), intent(in) :: coeff(:, :, :)
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(:)

        if (.not. is_setup) return
!$omp target exit data map(delete: y)
!$omp target exit data map(delete: x)
!$omp target exit data map(delete: coeff)
        is_setup = .false.
    end subroutine spline1d_many_omp_target_teardown

    subroutine spline1d_many_omp_target_eval_resident(order, num_points, &
                                                      num_quantities, &
                                                      periodic, x_min, h_step, &
                                                      coeff, x, y)
        integer, intent(in) :: order
        integer, intent(in) :: num_points
        integer, intent(in) :: num_quantities
        logical, intent(in) :: periodic
        real(dp), intent(in) :: x_min
        real(dp), intent(in) :: h_step
        real(dp), intent(in) :: coeff(num_quantities, 0:order, num_points)
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(num_quantities*size(x))

        integer :: ipt, iq, k_power, idx, base, k_wrap
        integer :: periodic_int
        real(dp) :: xj, x_norm, x_local, period, t, w

        if (periodic) then
            periodic_int = 1
        else
            periodic_int = 0
        end if
        period = h_step*real(num_points - 1, dp)

#include "spline1d_many_eval_body.inc"
    end subroutine spline1d_many_omp_target_eval_resident

end module spline1d_many_omp_target
