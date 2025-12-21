module spline1d_many_offload
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: spline1d_many_eval_host
    public :: spline1d_many_eval_resident
    public :: spline1d_many_teardown
    public :: spline1d_many_setup

    logical :: is_setup = .false.
contains

    subroutine spline1d_many_setup(coeff, x, y)
        real(dp), intent(in) :: coeff(:, :, :)
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(:)

        if (is_setup) return
#if defined(LIBNEO_ENABLE_OPENACC)
        !$acc enter data copyin(x)
        !$acc enter data create(y)
#endif
        is_setup = .true.
    end subroutine spline1d_many_setup

    subroutine spline1d_many_teardown(coeff, x, y)
        real(dp), intent(in) :: coeff(:, :, :)
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(:)

        if (.not. is_setup) return
#if defined(LIBNEO_ENABLE_OPENACC)
        !$acc exit data delete(y)
        !$acc exit data delete(x)
#endif
        is_setup = .false.
    end subroutine spline1d_many_teardown

    pure subroutine spline1d_many_eval_host(order, num_points, num_quantities, &
                                            periodic, x_min, h_step, coeff, x, y)
        integer, intent(in) :: order
        integer, intent(in) :: num_points
        integer, intent(in) :: num_quantities
        logical, intent(in) :: periodic
        real(dp), intent(in) :: x_min
        real(dp), intent(in) :: h_step
        real(dp), intent(in) :: coeff(num_quantities, 0:order, num_points)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(num_quantities*size(x))

        integer :: ipt, iq, k_power, idx, base, k_wrap
        integer :: periodic_int
        real(dp) :: xj, x_norm, x_local, period, t, w

        periodic_int = 0
        if (periodic) periodic_int = 1
        period = h_step*real(num_points - 1, dp)

        do ipt = 1, size(x)
            include "spline1d_many_point_body.inc"
        end do
    end subroutine spline1d_many_eval_host

    subroutine spline1d_many_eval_resident(order, num_points, num_quantities, &
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

        integer :: ipt, iq, k_power, idx, base, k_wrap
        integer :: periodic_int
        real(dp) :: xj, x_norm, x_local, period, t, w

        periodic_int = 0
        if (periodic) periodic_int = 1
        period = h_step*real(num_points - 1, dp)

#if defined(LIBNEO_ENABLE_OPENACC)
        !$acc parallel loop present(coeff, x, y) &
        !$acc& private(ipt, iq, k_power, idx, base, xj, x_norm, x_local, t, w, k_wrap) &
        !$acc& gang vector vector_length(256)
        do ipt = 1, size(x)
            include "spline1d_many_point_body.inc"
        end do
        !$acc end parallel loop
#else
        do ipt = 1, size(x)
            include "spline1d_many_point_body.inc"
        end do
#endif
    end subroutine spline1d_many_eval_resident

end module spline1d_many_offload
