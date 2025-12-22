module spline2d_many_offload
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: spline2d_many_setup
    public :: spline2d_many_teardown
    public :: spline2d_many_eval_host
    public :: spline2d_many_eval_resident

    logical :: is_setup = .false.
contains

    subroutine spline2d_many_setup(coeff, x, y)
        real(dp), intent(in) :: coeff(:, :, :, :, :)
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(inout) :: y(:)

        if (is_setup) return
#if defined(LIBNEO_ENABLE_OPENACC)
        !$acc enter data copyin(coeff, x)
        !$acc enter data create(y)
#endif
        is_setup = .true.
    end subroutine spline2d_many_setup

    subroutine spline2d_many_teardown(coeff, x, y)
        real(dp), intent(in) :: coeff(:, :, :, :, :)
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(inout) :: y(:)

        if (.not. is_setup) return
#if defined(LIBNEO_ENABLE_OPENACC)
        !$acc exit data delete(y)
        !$acc exit data delete(x)
        !$acc exit data delete(coeff)
#endif
        is_setup = .false.
    end subroutine spline2d_many_teardown

    pure subroutine spline2d_many_eval_host(order, num_points, num_quantities, &
                                            periodic, x_min, h_step, coeff, x, y)
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
        real(dp), intent(out) :: y(num_quantities*size(x, 2))

        integer :: ipt, iq, k1, k2, i1, i2, base, k_wrap
        integer :: periodic_int(2)
        integer :: npts, order1, order2
        real(dp) :: xj1, xj2, x_norm1, x_norm2, x_local1, x_local2
        real(dp) :: period(2)
        real(dp) :: t, w
        real(dp) :: v, yq

        if (size(x, 1) /= 2) error stop
        npts = size(x, 2)
        order1 = order(1)
        order2 = order(2)

        periodic_int = 0
        if (periodic(1)) periodic_int(1) = 1
        if (periodic(2)) periodic_int(2) = 1

        period(1) = h_step(1)*real(num_points(1) - 1, dp)
        period(2) = h_step(2)*real(num_points(2) - 1, dp)

        do ipt = 1, npts
            include "spline2d_many_point_body.inc"
        end do
    end subroutine spline2d_many_eval_host

    subroutine spline2d_many_eval_resident(order, num_points, num_quantities, &
                                           periodic, x_min, h_step, coeff, x, y)
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

        periodic_int = 0
        if (periodic(1)) periodic_int(1) = 1
        if (periodic(2)) periodic_int(2) = 1

        period(1) = h_step(1)*real(num_points(1) - 1, dp)
        period(2) = h_step(2)*real(num_points(2) - 1, dp)

#if defined(LIBNEO_ENABLE_OPENACC)
        !$acc parallel loop present(coeff, x, y) &
        !$acc& private(ipt, iq, k1, k2, i1, i2, base, xj1, xj2, x_norm1, x_norm2, &
        !$acc& x_local1, x_local2, t, w, v, yq, k_wrap) gang vector vector_length(256)
        do ipt = 1, npts
            include "spline2d_many_point_body.inc"
        end do
        !$acc end parallel loop
#else
        do ipt = 1, npts
            include "spline2d_many_point_body.inc"
        end do
#endif
    end subroutine spline2d_many_eval_resident

end module spline2d_many_offload
