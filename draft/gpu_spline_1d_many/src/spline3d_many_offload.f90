module spline3d_many_offload
    use, intrinsic :: iso_fortran_env, only: dp => real64
#if defined(LIBNEO_ENABLE_OPENACC)
    use openacc, only: acc_is_present
#endif
    implicit none
    private

    public :: spline3d_many_setup
    public :: spline3d_many_teardown
    public :: spline3d_many_eval_host
    public :: spline3d_many_eval_resident

    logical :: is_setup = .false.
#if defined(LIBNEO_ENABLE_OPENACC)
    logical :: did_copy_coeff = .false.
#endif

contains

    subroutine spline3d_many_setup(coeff, x, y)
        real(dp), intent(in) :: coeff(:, :, :, :, :, :, :)
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(inout) :: y(:)

        if (is_setup) return
#if defined(LIBNEO_ENABLE_OPENACC)
        did_copy_coeff = .false.
        if (.not. acc_is_present(coeff)) then
            !$acc enter data copyin(coeff)
            did_copy_coeff = .true.
        end if
        !$acc enter data copyin(x)
        !$acc enter data create(y)
#elif defined(LIBNEO_ENABLE_OPENMP)
        !$omp target enter data map(to: coeff, x) map(alloc: y)
#endif
        is_setup = .true.
    end subroutine spline3d_many_setup

    subroutine spline3d_many_teardown(coeff, x, y)
        real(dp), intent(in) :: coeff(:, :, :, :, :, :, :)
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(inout) :: y(:)

        if (.not. is_setup) return
#if defined(LIBNEO_ENABLE_OPENACC)
        !$acc exit data delete(y)
        !$acc exit data delete(x)
        if (did_copy_coeff) then
            !$acc exit data delete(coeff)
        end if
        did_copy_coeff = .false.
#elif defined(LIBNEO_ENABLE_OPENMP)
        !$omp target exit data map(delete: y, x, coeff)
#endif
        is_setup = .false.
    end subroutine spline3d_many_teardown

    pure subroutine spline3d_many_eval_host(order, num_points, num_quantities, &
                                            periodic, x_min, h_step, coeff, x, y)
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

        if (size(x, 1) /= 3) error stop
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

        do ipt = 1, npts
            include "spline3d_many_point_body.inc"
        end do
    end subroutine spline3d_many_eval_host

    subroutine spline3d_many_eval_resident(order, num_points, num_quantities, &
                                           periodic, x_min, h_step, coeff, x, y)
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
        real(dp), intent(inout) :: y(num_quantities*size(x, 2))

        integer, parameter :: threads_per_team = 256
        integer :: ipt, iq, k1, k2, k3, i1, i2, i3, base, k_wrap
        integer :: nteams
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

#if defined(LIBNEO_ENABLE_OPENACC)
        !$acc parallel loop present(coeff, x, y) &
        !$acc& private(ipt, iq, k1, k2, k3, i1, i2, i3, base, xj1, xj2, xj3) &
        !$acc& private(x_norm1, x_norm2, x_norm3, x_local1, x_local2, x_local3) &
        !$acc& private(t, w, v, w2, yq, k_wrap) gang vector vector_length(256)
        do ipt = 1, npts
            include "spline3d_many_point_body.inc"
        end do
        !$acc end parallel loop
#elif defined(LIBNEO_ENABLE_OPENMP)
        nteams = (npts + threads_per_team - 1)/threads_per_team
        !$omp target teams thread_limit(threads_per_team) num_teams(nteams)
        !$omp distribute parallel do simd &
        !$omp& private(iq, k1, k2, k3, i1, i2, i3, base, xj1, xj2, xj3) &
        !$omp& private(x_norm1, x_norm2, x_norm3, x_local1, x_local2, x_local3) &
        !$omp& private(t, w, v, w2, yq, k_wrap)
        do ipt = 1, npts
            include "spline3d_many_point_body.inc"
        end do
        !$omp end distribute parallel do simd
        !$omp end target teams
#else
        do ipt = 1, npts
            include "spline3d_many_point_body.inc"
        end do
#endif
    end subroutine spline3d_many_eval_resident

end module spline3d_many_offload
