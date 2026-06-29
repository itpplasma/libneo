!> Validate evaluate_batch_splines_3d_der3: its value, first and second
!> derivatives must match evaluate_batch_splines_3d_der2 to machine precision
!> (same spline), and its third derivatives must match central differences of
!> the second derivatives. Self-contained: builds a quintic spline of a smooth
!> analytic function, no external data.
program test_batch_interpolate_der3
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use batch_interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
        destroy_batch_splines_3d, evaluate_batch_splines_3d_der2, &
        evaluate_batch_splines_3d_der3
    implicit none

    integer, parameter :: NQ = 2, N1 = 24, N2 = 22, N3 = 20
    real(dp), parameter :: tol_match = 1.0e-9_dp, reltol_fd = 1.0e-4_dp, floor = 1.0e-6_dp
    real(dp) :: x_min(3), x_max(3), y(N1, N2, N3, NQ)
    integer :: order(3)
    logical :: periodic(3)
    type(BatchSplineData3D) :: spl
    real(dp) :: hgrid(3)
    integer :: i, j, k, q
    logical :: failed

    x_min = [0.0_dp, 0.0_dp, 0.0_dp]
    x_max = [1.0_dp, 1.0_dp, 1.0_dp]
    order = [5, 5, 5]
    periodic = [.false., .false., .false.]
    hgrid(1) = (x_max(1) - x_min(1))/real(N1 - 1, dp)
    hgrid(2) = (x_max(2) - x_min(2))/real(N2 - 1, dp)
    hgrid(3) = (x_max(3) - x_min(3))/real(N3 - 1, dp)

    do q = 1, NQ
        do k = 1, N3
            do j = 1, N2
                do i = 1, N1
                    y(i, j, k, q) = f(x_min(1) + real(i - 1, dp)*hgrid(1), &
                                      x_min(2) + real(j - 1, dp)*hgrid(2), &
                                      x_min(3) + real(k - 1, dp)*hgrid(3), q)
                end do
            end do
        end do
    end do

    call construct_batch_splines_3d(x_min, x_max, y, order, periodic, spl)

    failed = .false.
    call check_at([0.31_dp, 0.42_dp, 0.55_dp], failed)
    call check_at([0.62_dp, 0.18_dp, 0.73_dp], failed)
    call check_at([0.47_dp, 0.81_dp, 0.29_dp], failed)

    call destroy_batch_splines_3d(spl)

    if (failed) error stop "evaluate_batch_splines_3d_der3 mismatch"
    print *, "evaluate_batch_splines_3d_der3 matches der2 and finite differences"

contains

    real(dp) function f(a, b, c, q)
        real(dp), intent(in) :: a, b, c
        integer, intent(in) :: q
        f = sin(2.1_dp*a + 0.3_dp*q)*cos(1.7_dp*b)*exp(-0.4_dp*c) &
            + 0.5_dp*a*a*b*c
    end function f

    subroutine eval_d2(x, d2)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: d2(6, NQ)
        real(dp) :: yv(NQ), dy(3, NQ)
        call evaluate_batch_splines_3d_der2(spl, x, yv, dy, d2)
    end subroutine eval_d2

    subroutine check_at(x, failed)
        real(dp), intent(in) :: x(3)
        logical, intent(inout) :: failed
        real(dp) :: yv(NQ), dy(3, NQ), d2(6, NQ), d3(10, NQ)
        real(dp) :: y2(NQ), dy2(3, NQ), d2b(6, NQ)
        real(dp) :: fd1(6), fd2(6), fd3(6), fd(10, NQ)
        real(dp), parameter :: h = 1.0e-3_dp
        integer :: q, m

        call evaluate_batch_splines_3d_der3(spl, x, yv, dy, d2, d3)
        call evaluate_batch_splines_3d_der2(spl, x, y2, dy2, d2b)

        ! value, first, second derivatives must match der2 to machine precision.
        do q = 1, NQ
            call near("val", abs(yv(q) - y2(q)), failed)
            do m = 1, 3
                call near("d1", abs(dy(m, q) - dy2(m, q)), failed)
            end do
            do m = 1, 6
                call near("d2", abs(d2(m, q) - d2b(m, q)), failed)
            end do
        end do

        ! third derivatives vs central differences of d2 in x1, x2, x3.
        do q = 1, NQ
            call d2deriv(x, h, 1, q, fd1)
            call d2deriv(x, h, 2, q, fd2)
            call d2deriv(x, h, 3, q, fd3)
            fd(1, q) = fd1(1)   ! 111 = d(11)/d1
            fd(2, q) = fd2(1)   ! 112 = d(11)/d2
            fd(3, q) = fd3(1)   ! 113 = d(11)/d3
            fd(4, q) = fd1(4)   ! 122 = d(22)/d1
            fd(5, q) = fd3(2)   ! 123 = d(12)/d3
            fd(6, q) = fd1(6)   ! 133 = d(33)/d1
            fd(7, q) = fd2(4)   ! 222 = d(22)/d2
            fd(8, q) = fd3(4)   ! 223 = d(22)/d3
            fd(9, q) = fd2(6)   ! 233 = d(33)/d2
            fd(10, q) = fd3(6)  ! 333 = d(33)/d3
            do m = 1, 10
                call near_rel("d3", d3(m, q), fd(m, q), failed)
            end do
        end do
    end subroutine check_at

    subroutine d2deriv(x, h, dir, q, out)
        real(dp), intent(in) :: x(3), h
        integer, intent(in) :: dir, q
        real(dp), intent(out) :: out(6)
        real(dp) :: xp(3), dp2(6, NQ), dm2(6, NQ)
        xp = x; xp(dir) = x(dir) + h; call eval_d2(xp, dp2)
        xp = x; xp(dir) = x(dir) - h; call eval_d2(xp, dm2)
        out = (dp2(:, q) - dm2(:, q))/(2.0_dp*h)
    end subroutine d2deriv

    subroutine near(name, diff, failed)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: diff
        logical, intent(inout) :: failed
        if (diff > tol_match) then
            print '(a,a,1x,es12.4)', name, " der3-vs-der2 mismatch ", diff
            failed = .true.
        end if
    end subroutine near

    subroutine near_rel(name, a, b, failed)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a, b
        logical, intent(inout) :: failed
        real(dp) :: denom
        denom = max(abs(a), abs(b), floor)
        if (abs(a - b) > reltol_fd*denom) then
            print '(a,1x,a,es15.7,1x,a,es15.7,1x,es11.3)', name, "analytic=", a, &
                "fd=", b, abs(a - b)/denom
            failed = .true.
        end if
    end subroutine near_rel

end program test_batch_interpolate_der3
