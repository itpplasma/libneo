module spline1d_many_openacc
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: spline1d_many_openacc_eval_resident
    public :: spline1d_many_openacc_teardown
    public :: spline1d_many_openacc_setup

    logical :: is_setup = .false.

contains

    subroutine spline1d_many_openacc_setup(coeff, x, y)
        real(dp), intent(in) :: coeff(:, :, :)
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(:)

        if (is_setup) return
        !$acc enter data copyin(coeff)
        !$acc enter data copyin(x)
        !$acc enter data create(y)
        is_setup = .true.
    end subroutine spline1d_many_openacc_setup

    subroutine spline1d_many_openacc_teardown(coeff, x, y)
        real(dp), intent(in) :: coeff(:, :, :)
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(:)

        if (.not. is_setup) return
        !$acc exit data delete(y)
        !$acc exit data delete(x)
        !$acc exit data delete(coeff)
        is_setup = .false.
    end subroutine spline1d_many_openacc_teardown

    subroutine spline1d_many_openacc_eval_resident(order, num_points, num_quantities, &
                                                   periodic, x_min, h_step, coeff, x, y)
        integer, intent(in) :: order
        integer, intent(in) :: num_points
        integer, intent(in) :: num_quantities
        logical, intent(in) :: periodic
        real(dp), intent(in) :: x_min
        real(dp), intent(in) :: h_step
        real(dp), intent(in) :: coeff(num_quantities, 0:order, num_points)
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(num_quantities*size(x))

        integer :: ipt, iq, k_power, idx, base
        real(dp) :: xj, x_norm, x_local, period
        logical :: is_periodic

        is_periodic = periodic
        period = h_step*real(num_points - 1, dp)

        !$acc parallel loop present(coeff, x, y) &
        !$acc& private(xj, x_norm, x_local, idx, base) gang
        do ipt = 1, size(x)
            base = (ipt - 1)*num_quantities
            if (is_periodic) then
                xj = modulo(x(ipt) - x_min, period) + x_min
            else
                xj = x(ipt)
            end if

            x_norm = (xj - x_min)/h_step
            idx = max(0, min(num_points - 2, int(x_norm)))
            x_local = (x_norm - real(idx, dp))*h_step

            !$acc loop vector
            do iq = 1, num_quantities
                y(base + iq) = coeff(iq, order, idx + 1)
            end do

            do k_power = order - 1, 0, -1
                !$acc loop vector
                do iq = 1, num_quantities
                    y(base + iq) = coeff(iq, k_power, idx + 1) + x_local*y(base + iq)
                end do
            end do
        end do
        !$acc end parallel loop
    end subroutine spline1d_many_openacc_eval_resident

end module spline1d_many_openacc
