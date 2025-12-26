program test_batch_interpolate_mask
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use batch_interpolate, only: BatchSplineData1D, BatchSplineData3D
    use batch_interpolate, only: construct_batch_splines_1d, destroy_batch_splines_1d
    use batch_interpolate, only: construct_batch_splines_3d, destroy_batch_splines_3d
    use batch_interpolate, only: evaluate_batch_splines_1d_many_der3
    use batch_interpolate, only: evaluate_batch_splines_3d_many_der2
    use batch_interpolate, only: evaluate_batch_splines_1d_many_der3_mask
    use batch_interpolate, only: evaluate_batch_splines_3d_many_der2_mask
    implicit none

    call test_1d_der3_mask()
    call test_3d_der2_mask()

contains

    subroutine assert_close(name, a, b, tol)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a(:)
        real(dp), intent(in) :: b(:)
        real(dp), intent(in) :: tol
        real(dp) :: diff

        diff = maxval(abs(a - b))
        if (diff > tol) then
            write (*, *) "FAIL", trim(name), "diff", diff, "tol", tol
            error stop 1
        end if
    end subroutine assert_close


    subroutine assert_unchanged(name, a, sentinel)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a(:)
        real(dp), intent(in) :: sentinel

        if (any(a /= sentinel)) then
            write (*, *) "FAIL", trim(name), "unexpected overwrite"
            error stop 1
        end if
    end subroutine assert_unchanged


    subroutine test_1d_der3_mask()
        integer, parameter :: order = 5
        integer, parameter :: ngrid = 32
        integer, parameter :: nq = 3
        integer, parameter :: npts = 17
        logical, parameter :: periodic = .true.
        real(dp), parameter :: x_min = 0.0d0
        real(dp), parameter :: x_max = 2.0d0
        real(dp), parameter :: tol = 1.0d-12
        real(dp), parameter :: sentinel = -999.0d0

        type(BatchSplineData1D) :: spl
        real(dp) :: x_grid(ngrid)
        real(dp) :: y_grid(ngrid, nq)
        real(dp) :: x_eval(npts)
        logical :: mask(npts)

        real(dp) :: y_full(nq, npts), dy_full(nq, npts), d2_full(nq, npts), d3_full(nq, npts)
        real(dp) :: y_mask(nq, npts), dy_mask(nq, npts), d2_mask(nq, npts), d3_mask(nq, npts)
        integer :: i, iq

        do i = 1, ngrid
            x_grid(i) = x_min + (x_max - x_min) * real(i - 1, dp) / real(ngrid - 1, dp)
        end do
        do iq = 1, nq
            do i = 1, ngrid
                y_grid(i, iq) = cos(x_grid(i) + 0.1d0 * real(iq - 1, dp))
            end do
        end do
        call construct_batch_splines_1d(x_min, x_max, y_grid, order, periodic, spl)

        do i = 1, npts
            x_eval(i) = x_min + 0.73d0 * real(i - 1, dp)
            mask(i) = mod(i, 2) == 1
        end do

        call evaluate_batch_splines_1d_many_der3(spl, x_eval, y_full, dy_full, d2_full, d3_full)

        y_mask = sentinel
        dy_mask = sentinel
        d2_mask = sentinel
        d3_mask = sentinel
        call evaluate_batch_splines_1d_many_der3_mask(spl, x_eval, mask, y_mask, dy_mask, &
                                                      d2_mask, d3_mask)

        do i = 1, npts
            if (mask(i)) then
                call assert_close("1d_y", y_full(:, i), y_mask(:, i), tol)
                call assert_close("1d_dy", dy_full(:, i), dy_mask(:, i), tol)
                call assert_close("1d_d2", d2_full(:, i), d2_mask(:, i), tol)
                call assert_close("1d_d3", d3_full(:, i), d3_mask(:, i), tol)
            else
                call assert_unchanged("1d_y_sentinel", y_mask(:, i), sentinel)
                call assert_unchanged("1d_dy_sentinel", dy_mask(:, i), sentinel)
                call assert_unchanged("1d_d2_sentinel", d2_mask(:, i), sentinel)
                call assert_unchanged("1d_d3_sentinel", d3_mask(:, i), sentinel)
            end if
        end do

        call destroy_batch_splines_1d(spl)
    end subroutine test_1d_der3_mask


    subroutine test_3d_der2_mask()
        integer, parameter :: order(3) = [5, 5, 5]
        integer, parameter :: ngrid(3) = [16, 12, 10]
        integer, parameter :: nq = 2
        integer, parameter :: npts = 19
        logical, parameter :: periodic(3) = [.false., .true., .true.]
        real(dp), parameter :: x_min(3) = [0.0d0, 0.0d0, 0.0d0]
        real(dp), parameter :: x_max(3) = [1.0d0, 2.0d0, 3.0d0]
        real(dp), parameter :: tol = 1.0d-12
        real(dp), parameter :: sentinel = -999.0d0

        type(BatchSplineData3D) :: spl
        real(dp) :: y_grid(ngrid(1), ngrid(2), ngrid(3), nq)
        real(dp) :: x_eval(3, npts)
        logical :: mask(npts)

        real(dp) :: y_full(nq, npts), y_mask(nq, npts)
        real(dp) :: dy_full(3, nq, npts), dy_mask(3, nq, npts)
        real(dp) :: d2_full(6, nq, npts), d2_mask(6, nq, npts)
        integer :: i1, i2, i3, iq, ip

        do iq = 1, nq
            do i3 = 1, ngrid(3)
                do i2 = 1, ngrid(2)
                    do i1 = 1, ngrid(1)
                        y_grid(i1, i2, i3, iq) = &
                            cos(x_min(1) + real(i1 - 1, dp) * (x_max(1) - x_min(1)) / &
                                real(ngrid(1) - 1, dp)) * &
                            cos(x_min(2) + real(i2 - 1, dp) * (x_max(2) - x_min(2)) / &
                                real(ngrid(2) - 1, dp)) * &
                            cos(x_min(3) + real(i3 - 1, dp) * (x_max(3) - x_min(3)) / &
                                real(ngrid(3) - 1, dp)) * &
                            (1.0d0 + 0.07d0 * real(iq - 1, dp))
                    end do
                end do
            end do
        end do

        call construct_batch_splines_3d(x_min, x_max, y_grid, order, periodic, spl)

        do ip = 1, npts
            x_eval(1, ip) = x_min(1) + (x_max(1) - x_min(1)) * real(ip - 1, dp) / real(npts, dp)
            x_eval(2, ip) = x_min(2) + (x_max(2) - x_min(2)) * real(2 * ip - 1, dp) / &
                            real(2 * npts, dp)
            x_eval(3, ip) = x_min(3) + (x_max(3) - x_min(3)) * real(3 * ip - 1, dp) / &
                            real(3 * npts, dp)
            mask(ip) = mod(ip, 3) /= 0
        end do

        call evaluate_batch_splines_3d_many_der2(spl, x_eval, y_full, dy_full, d2_full)

        y_mask = sentinel
        dy_mask = sentinel
        d2_mask = sentinel
        call evaluate_batch_splines_3d_many_der2_mask(spl, x_eval, mask, y_mask, dy_mask, d2_mask)

        do ip = 1, npts
            if (mask(ip)) then
                call assert_close("3d_y", y_full(:, ip), y_mask(:, ip), tol)
                call assert_close("3d_dy", reshape(dy_full(:, :, ip), [3 * nq]), &
                                  reshape(dy_mask(:, :, ip), [3 * nq]), tol)
                call assert_close("3d_d2", reshape(d2_full(:, :, ip), [6 * nq]), &
                                  reshape(d2_mask(:, :, ip), [6 * nq]), tol)
            else
                call assert_unchanged("3d_y_sentinel", y_mask(:, ip), sentinel)
                call assert_unchanged("3d_dy_sentinel", reshape(dy_mask(:, :, ip), [3 * nq]), &
                                      sentinel)
                call assert_unchanged("3d_d2_sentinel", reshape(d2_mask(:, :, ip), [6 * nq]), &
                                      sentinel)
            end if
        end do

        call destroy_batch_splines_3d(spl)
    end subroutine test_3d_der2_mask

end program test_batch_interpolate_mask
