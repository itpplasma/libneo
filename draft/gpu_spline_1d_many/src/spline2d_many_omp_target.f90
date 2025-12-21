module spline2d_many_omp_target
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: spline2d_many_omp_target_eval_resident

contains

    subroutine spline2d_many_omp_target_eval_resident(order, num_points, &
                                                      num_quantities, &
                                                      periodic, x_min, h_step, &
                                                      coeff, x, y)
        integer, intent(in) :: order(2)
        integer, intent(in) :: num_points(2)
        integer, intent(in) :: num_quantities
        logical, intent(in) :: periodic(2)
        real(dp), intent(in) :: x_min(2)
        real(dp), intent(in) :: h_step(2)
        real(dp), intent(in) :: coeff(num_quantities, 0:order(1), 0:order(2), &
                                      num_points(1), &
                                      num_points(2))
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(inout) :: y(num_quantities*size(x, 2))

        integer :: ipt, iq, k1, k2, i1, i2, base, k_wrap
        integer :: periodic_int(2)
        integer :: npts, order1, order2
        real(dp) :: xj1, xj2, x_norm1, x_norm2, x_local1, x_local2
        real(dp) :: period(2)
        real(dp) :: t, w
        real(dp) :: v, yq

        npts = size(x, 2)
        order1 = order(1)
        order2 = order(2)

        if (periodic(1)) then
            periodic_int(1) = 1
        else
            periodic_int(1) = 0
        end if
        if (periodic(2)) then
            periodic_int(2) = 1
        else
            periodic_int(2) = 0
        end if

        period(1) = h_step(1)*real(num_points(1) - 1, dp)
        period(2) = h_step(2)*real(num_points(2) - 1, dp)

#include "spline2d_many_eval_body.inc"
    end subroutine spline2d_many_omp_target_eval_resident

end module spline2d_many_omp_target
