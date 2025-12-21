module spline1d_many_cpu
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: spline1d_many_cpu_eval

contains

    pure subroutine spline1d_many_cpu_eval(order, num_points, num_quantities, &
                                           periodic, &
                                           x_min, h_step, coeff, x, y)
        integer, intent(in) :: order
        integer, intent(in) :: num_points
        integer, intent(in) :: num_quantities
        logical, intent(in) :: periodic
        real(dp), intent(in) :: x_min
        real(dp), intent(in) :: h_step
        real(dp), intent(in) :: coeff(num_quantities, 0:order, num_points)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(num_quantities*size(x))

        integer :: ipt, iq, k_power, idx, base
        integer :: periodic_int
        real(dp) :: xj, x_norm, x_local, period, t, w

        period = h_step*real(num_points - 1, dp)
        if (periodic) then
            periodic_int = 1
        else
            periodic_int = 0
        end if

#include "spline1d_many_eval_body.inc"
    end subroutine spline1d_many_cpu_eval

end module spline1d_many_cpu
