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

        integer :: ipt, iq, k_power, idx
        real(dp) :: xj, x_norm, x_local, period
        integer :: base

        period = h_step*real(num_points - 1, dp)

        do ipt = 1, size(x)
            if (periodic) then
                xj = modulo(x(ipt) - x_min, period) + x_min
            else
                xj = x(ipt)
            end if

            x_norm = (xj - x_min)/h_step
            idx = max(0, min(num_points - 2, int(x_norm)))
            x_local = (x_norm - real(idx, dp))*h_step

            base = (ipt - 1)*num_quantities
            do iq = 1, num_quantities
                y(base + iq) = coeff(iq, order, idx + 1)
            end do

            do k_power = order - 1, 0, -1
                do iq = 1, num_quantities
                    y(base + iq) = coeff(iq, k_power, idx + 1) + x_local*y(base + iq)
                end do
            end do
        end do
    end subroutine spline1d_many_cpu_eval

end module spline1d_many_cpu
