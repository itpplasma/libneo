module spline3d_many_cpu
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: spline3d_many_cpu_eval

contains

    pure subroutine spline3d_many_cpu_eval(order, num_points, num_quantities, &
                                           periodic, x_min, &
                                           h_step, coeff, x, y)
        integer, intent(in) :: order(3)
        integer, intent(in) :: num_points(3)
        integer, intent(in) :: num_quantities
        logical, intent(in) :: periodic(3)
        real(dp), intent(in) :: x_min(3)
        real(dp), intent(in) :: h_step(3)
        real(dp), intent(in) :: coeff(num_quantities, 0:order(1), 0:order(2), &
                                      0:order(3), &
                                      num_points(1), num_points(2), num_points(3))
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: y(num_quantities*size(x, 2))

        integer :: ipt, iq, k1, k2, k3, i1, i2, i3, base, k_wrap
        integer :: periodic_int(3)
        integer :: npts, order1, order2, order3
        real(dp) :: xj1, xj2, xj3, x_norm1, x_norm2, x_norm3
        real(dp) :: x_local1, x_local2, x_local3
        real(dp) :: period(3)
        real(dp) :: t, w
        real(dp) :: v, w2, yq

        npts = size(x, 2)
        order1 = order(1)
        order2 = order(2)
        order3 = order(3)

        periodic_int = 0
        if (periodic(1)) periodic_int(1) = 1
        if (periodic(2)) periodic_int(2) = 1
        if (periodic(3)) periodic_int(3) = 1

        period(1) = h_step(1)*real(num_points(1) - 1, dp)
        period(2) = h_step(2)*real(num_points(2) - 1, dp)
        period(3) = h_step(3)*real(num_points(3) - 1, dp)

#include "spline3d_many_eval_body.inc"
    end subroutine spline3d_many_cpu_eval

end module spline3d_many_cpu

