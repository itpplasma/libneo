program test_batch_interpolate_o333
    ! Guards the cubic path of the 3D batch-spline evaluators (generated cell
    ! kernels for value, first and second derivatives).
    !
    ! The specialisation exists for speed: chartmap coordinates are cubic
    ! splines of the position, evaluated with second derivatives at every
    ! orbit step of SIMPLE's reference-coordinate tracing. The oracle is the
    ! general kernel of each family, which the generated kernels must match
    ! to 1e-12 relative (they order the Horner recurrences differently, so
    ! not bitwise). Cases cover 1, 3 (chart position) and MAX_QUANTITIES
    ! quantities, all-periodic axes and a non-periodic radial axis.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use batch_interpolate, only: BatchSplineData3D, construct_batch_splines_3d, &
                                 destroy_batch_splines_3d, &
                                 evaluate_batch_splines_3d, &
                                 evaluate_batch_splines_3d_many, &
                                 evaluate_batch_splines_3d_der, &
                                 evaluate_batch_splines_3d_der_general, &
                                 evaluate_batch_splines_3d_der2, &
                                 evaluate_batch_splines_3d_der2_general
    implicit none

    integer :: nfail
    nfail = 0
    call check_case(1, [.false., .true., .true.], nfail)
    call check_case(3, [.false., .true., .true.], nfail)
    call check_case(3, [.true., .true., .true.], nfail)
    call check_case(8, [.false., .true., .true.], nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "test_batch_interpolate_o333: all tests passed"

contains

    real(dp) function sample_field(x1, x2, x3, iq) result(v)
        real(dp), intent(in) :: x1, x2, x3
        integer, intent(in) :: iq
        v = cos(x1 + 0.1_dp*iq)*cos(x2)*cos(x3) + 0.25_dp*sin(2.0_dp*x1)*cos(x3) &
            + 0.1_dp*real(iq - 1, dp)*sin(x2 + x3)
    end function sample_field

    real(dp) function rel(a, b)
        real(dp), intent(in) :: a(:), b(:)
        rel = maxval(abs(a - b)/max(1.0_dp, abs(b)))
    end function rel

    subroutine check_case(nq, periodic, nfail)
        integer, intent(in) :: nq
        logical, intent(in) :: periodic(3)
        integer, intent(inout) :: nfail

        integer, parameter :: order(3) = [3, 3, 3]
        integer, parameter :: ngrid(3) = [14, 12, 11]
        integer, parameter :: npts = 500
        real(dp), parameter :: x_min(3) = [0.0_dp, 0.0_dp, 0.0_dp]
        real(dp), parameter :: x_max(3) = [1.7_dp, 2.3_dp, 1.1_dp]
        real(dp), parameter :: tol = 1.0e-12_dp

        type(BatchSplineData3D) :: spl
        real(dp), allocatable :: y(:, :, :, :)
        real(dp) :: x(3), p(3), u, w0, w1, w2
        real(dp) :: v(nq), vr(nq, 1), dv(3, nq), dvr(3, nq), v1(nq), v1r(nq)
        real(dp) :: v2(nq), dv2(3, nq), d2v(6, nq), v2r(nq), dv2r(3, nq), d2vr(6, nq)
        integer :: i1, i2, i3, iq, ip, j

        allocate (y(ngrid(1), ngrid(2), ngrid(3), nq))
        do iq = 1, nq
            do i3 = 1, ngrid(3)
                do i2 = 1, ngrid(2)
                    do i1 = 1, ngrid(1)
                        p = x_min + real([i1, i2, i3] - 1, dp)*(x_max - x_min) &
                            /real(ngrid - 1, dp)
                        y(i1, i2, i3, iq) = sample_field(p(1), p(2), p(3), iq)
                    end do
                end do
            end do
        end do
        call construct_batch_splines_3d(x_min, x_max, y, order, periodic, spl)

        w0 = 0.0_dp; w1 = 0.0_dp; w2 = 0.0_dp
        do ip = 1, npts
            do j = 1, 3
                u = mod(real(ip, dp)*(0.6180339887498949_dp + 0.11_dp*real(j, dp)), &
                        1.0_dp)
                if (periodic(j)) then
                    x(j) = x_min(j) - (x_max(j) - x_min(j)) + &
                           3.0_dp*(x_max(j) - x_min(j))*u
                else
                    x(j) = x_min(j) + (x_max(j) - x_min(j))*u
                end if
            end do
            call evaluate_batch_splines_3d(spl, x, v)
            call evaluate_batch_splines_3d_many(spl, reshape(x, [3, 1]), vr)
            w0 = max(w0, rel(v, vr(:, 1)))

            call evaluate_batch_splines_3d_der(spl, x, v1, dv)
            call evaluate_batch_splines_3d_der_general(spl, x, v1r, dvr)
            w1 = max(w1, rel(v1, v1r), rel(reshape(dv, [3*nq]), reshape(dvr, [3*nq])))

            call evaluate_batch_splines_3d_der2(spl, x, v2, dv2, d2v)
            call evaluate_batch_splines_3d_der2_general(spl, x, v2r, dv2r, d2vr)
            w2 = max(w2, rel(v2, v2r), rel(reshape(dv2, [3*nq]), reshape(dv2r, [3*nq])), &
                     rel(reshape(d2v, [6*nq]), reshape(d2vr, [6*nq])))
        end do

        write (*, "(a,i0,a,3l2,a,3es10.3)") "  nq=", nq, " periodic=", periodic, &
            ": worst rel diff to the general kernels (value, d1, d2) ", w0, w1, w2
        if (.not. (w0 < tol .and. w1 < tol .and. w2 < tol)) then
            write (error_unit, "(a)") "  cubic kernels differ from the general kernels"
            nfail = nfail + 1
        end if

        call destroy_batch_splines_3d(spl)
    end subroutine check_case

end program test_batch_interpolate_o333
