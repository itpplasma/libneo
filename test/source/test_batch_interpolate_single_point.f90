program test_batch_interpolate_single_point
    ! Guards the single-point fast paths used in orbit tracing.
    !
    !   1. evaluate_batch_splines_1d_der2/der3 call their single-point cores
    !      directly. They must equal the many-point routines, which loop over
    !      the same cores: bitwise with the argument "bitwise"
    !      (LIBNEO_DETERMINISTIC_FP builds), else to 1e-13 relative.
    !   2. evaluate_batch_splines_3d_der2 sends single-quantity quintic splines
    !      to the nq1_o555 kernel. Value, first and second derivatives must
    !      agree with evaluate_batch_splines_3d_der3, a separate kernel family.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use batch_interpolate, only: BatchSplineData1D, BatchSplineData3D, &
                                 construct_batch_splines_1d, &
                                 destroy_batch_splines_1d, &
                                 evaluate_batch_splines_1d_der2, &
                                 evaluate_batch_splines_1d_der3, &
                                 evaluate_batch_splines_1d_many_der2, &
                                 evaluate_batch_splines_1d_many_der3, &
                                 construct_batch_splines_3d, &
                                 destroy_batch_splines_3d, &
                                 evaluate_batch_splines_3d_der2, &
                                 evaluate_batch_splines_3d_der3
    implicit none

    integer :: nfail
    logical :: bitwise
    character(len=16) :: arg

    call get_command_argument(1, arg)
    bitwise = trim(arg) == "bitwise"
    nfail = 0
    call check_1d(.false., nfail)
    call check_1d(.true., nfail)
    call check_3d_der2_nq1_o555(.false., nfail)
    call check_3d_der2_nq1_o555(.true., nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "test_batch_interpolate_single_point: all tests passed"

contains

    subroutine check_1d(periodic, nfail)
        logical, intent(in) :: periodic
        integer, intent(inout) :: nfail

        integer, parameter :: n = 23, nq = 3, npts = 300
        type(BatchSplineData1D) :: spl
        real(dp) :: y(n, nq), x(1), u, worst
        real(dp) :: v(nq), d1(nq), d2(nq), d3(nq)
        real(dp) :: vm(nq, 1), d1m(nq, 1), d2m(nq, 1), d3m(nq, 1)
        integer :: i, iq, ip, order, nbad

        do iq = 1, nq
            do i = 1, n
                u = 2.0_dp*real(i - 1, dp)/real(n - 1, dp)
                y(i, iq) = sin(3.0_dp*u + iq) + 0.2_dp*iq*cos(u)
            end do
        end do
        nbad = 0
        worst = 0.0_dp
        do order = 3, 5, 2
            call construct_batch_splines_1d(0.0_dp, 2.0_dp, y, order, periodic, spl)
            do ip = 1, npts
                x(1) = -0.5_dp + 3.0_dp*mod(real(ip, dp)*0.6180339887498949_dp, 1.0_dp)
                if (.not. periodic) x(1) = 2.0_dp*mod(real(ip, dp)*0.618034_dp, 1.0_dp)
                call evaluate_batch_splines_1d_der2(spl, x(1), v, d1, d2)
                call evaluate_batch_splines_1d_many_der2(spl, x, vm, d1m, d2m)
                nbad = nbad + count(v /= vm(:, 1)) + count(d1 /= d1m(:, 1)) &
                       + count(d2 /= d2m(:, 1))
                call evaluate_batch_splines_1d_der3(spl, x(1), v, d1, d2, d3)
                call evaluate_batch_splines_1d_many_der3(spl, x, vm, d1m, d2m, d3m)
                nbad = nbad + count(v /= vm(:, 1)) + count(d1 /= d1m(:, 1)) &
                       + count(d2 /= d2m(:, 1)) + count(d3 /= d3m(:, 1))
                worst = max(worst, rel(v, vm(:, 1)), rel(d1, d1m(:, 1)), &
                            rel(d2, d2m(:, 1)), rel(d3, d3m(:, 1)))
            end do
            call destroy_batch_splines_1d(spl)
        end do
        write (*, "(a,l2,a,i0,a,es10.3)") "  1d periodic=", periodic, &
            ": single- vs many-point mismatches ", nbad, ", worst rel ", worst
        if ((bitwise .and. nbad > 0) .or. .not. (worst < 1.0e-13_dp)) &
            nfail = nfail + 1
    end subroutine check_1d

    pure real(dp) function rel(a, b)
        real(dp), intent(in) :: a(:), b(:)
        rel = maxval(abs(a - b)/max(1.0_dp, abs(b)))
    end function rel

    subroutine check_3d_der2_nq1_o555(periodic1, nfail)
        logical, intent(in) :: periodic1
        integer, intent(inout) :: nfail

        integer, parameter :: ngrid(3) = [14, 12, 11], npts = 500
        real(dp), parameter :: x_min(3) = 0.0_dp, x_max(3) = [1.7_dp, 2.3_dp, 1.1_dp]
        real(dp), parameter :: tol = 1.0e-12_dp
        type(BatchSplineData3D) :: spl
        real(dp), allocatable :: y(:, :, :, :)
        real(dp) :: x(3), p(3), u, worst
        real(dp) :: v(1), d1(3, 1), d2(6, 1)
        real(dp) :: v3(1), d13(3, 1), d23(6, 1), d33(10, 1)
        integer :: i1, i2, i3, ip, j

        allocate (y(ngrid(1), ngrid(2), ngrid(3), 1))
        do i3 = 1, ngrid(3)
            do i2 = 1, ngrid(2)
                do i1 = 1, ngrid(1)
                    p = x_min + real([i1, i2, i3] - 1, dp)*(x_max - x_min) &
                        /real(ngrid - 1, dp)
                    y(i1, i2, i3, 1) = cos(p(1))*cos(p(2))*cos(p(3)) &
                                       + 0.25_dp*sin(2.0_dp*p(1))*cos(p(3))
                end do
            end do
        end do
        call construct_batch_splines_3d(x_min, x_max, y, [5, 5, 5], &
                                        [periodic1, .true., .true.], spl)
        worst = 0.0_dp
        do ip = 1, npts
            do j = 1, 3
                u = mod(real(ip, dp)*(0.6180339887498949_dp + 0.11_dp*j), 1.0_dp)
                x(j) = x_min(j) + (x_max(j) - x_min(j))*u
            end do
            call evaluate_batch_splines_3d_der2(spl, x, v, d1, d2)
            call evaluate_batch_splines_3d_der3(spl, x, v3, d13, d23, d33)
            worst = max(worst, abs(v(1) - v3(1)), &
                        maxval(abs(d1 - d13)/max(1.0_dp, abs(d13))), &
                        maxval(abs(d2 - d23)/max(1.0_dp, abs(d23))))
        end do
        write (*, "(a,l2,a,es10.3)") "  3d der2 nq1 o555 periodic1=", periodic1, &
            ": worst rel diff vs der3 ", worst
        if (.not. (worst < tol)) nfail = nfail + 1
        call destroy_batch_splines_3d(spl)
    end subroutine check_3d_der2_nq1_o555

end program test_batch_interpolate_single_point
