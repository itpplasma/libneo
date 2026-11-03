program test_batch_interpolate_horner
    ! Quintic batch splines evaluate through generated cell kernels. Every value
    ! and derivative they return is compared with a dense sum over the cell's
    ! monomial coefficients, evaluated here in extended precision:
    !   d^a p = sum_k c_k prod_j k_j!/(k_j - a_j)! x_j**(k_j - a_j),
    ! with the cell and local coordinate located independently on a dyadic
    ! grid. The bound is 64 ulp of the summed term magnitudes, so a wrong
    ! derivative term fails.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use batch_interpolate, only: BatchSplineData1D, BatchSplineData3D, &
                                 construct_batch_splines_1d, &
                                 destroy_batch_splines_1d, &
                                 evaluate_batch_splines_1d_der2, &
                                 evaluate_batch_splines_1d_der3, &
                                 evaluate_batch_spline_1d_scalar_quintic_der2, &
                                 evaluate_batch_spline_1d_pair_quintic_der2, &
                                 construct_batch_splines_3d, &
                                 destroy_batch_splines_3d, &
                                 evaluate_batch_splines_3d, &
                                 evaluate_batch_splines_3d_der, &
                                 evaluate_batch_splines_3d_der2, &
                                 evaluate_batch_splines_3d_der2_rmix, &
                                 evaluate_batch_spline_3d_scalar_quintic_der2_rmix
    implicit none

    integer, parameter :: qp = selected_real_kind(30)
    integer, parameter :: hessian(3, 6) = reshape([2, 0, 0, 1, 1, 0, 1, 0, 1, &
                                                   0, 2, 0, 0, 1, 1, 0, 0, 2], [3, 6])
    integer :: nfail, nq

    nfail = 0
    do nq = 1, 7
        call check_3d(nq, [.false., .true., .true.], nfail)
    end do
    call check_3d(3, [.true., .true., .true.], nfail)
    call check_1d(1, .false., nfail)
    call check_1d(2, .false., nfail)
    call check_1d(3, .true., nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " check(s) failed"
        stop 1
    end if
    write (*, "(a)") "test_batch_interpolate_horner: all tests passed"

contains

    subroutine check_3d(nq, periodic, nfail)
        integer, intent(in) :: nq
        logical, intent(in) :: periodic(3)
        integer, intent(inout) :: nfail

        integer, parameter :: ngrid(3) = [9, 8, 7], npts = 200
        real(dp), parameter :: x_min(3) = [0.0_dp, -1.0_dp, 0.5_dp]
        real(dp), parameter :: x_max(3) = x_min + real(ngrid - 1, dp)*0.125_dp
        type(BatchSplineData3D) :: spl
        real(dp), allocatable :: y(:, :, :, :)
        real(dp) :: x(3), u(3), xl(3)
        real(dp) :: v(nq), dv(3, nq), v2(nq), dv2(3, nq), d2v(6, nq)
        real(dp) :: v3(nq), dv3(3, nq), d2r(3, nq), v0(nq)
        real(dp) :: vs, dvs(3), d2s(3)
        integer :: ip, iq, j, cell(3), nbad

        call random_seed(put=[(31*j + nq, j=1, 64)])
        allocate (y(ngrid(1), ngrid(2), ngrid(3), nq))
        call random_number(y)
        call construct_batch_splines_3d(x_min, x_max, y, [5, 5, 5], periodic, spl)
        nbad = 0
        do ip = 1, npts
            call random_number(u)
            x = x_min + dyadic(u)*(x_max - x_min)
            if (ip <= 3) x = x_min + real(ip - 1, dp)*0.4375_dp*(x_max - x_min)
            call evaluate_batch_splines_3d(spl, x, v0)
            call evaluate_batch_splines_3d_der(spl, x, v, dv)
            call evaluate_batch_splines_3d_der2(spl, x, v2, dv2, d2v)
            call evaluate_batch_splines_3d_der2_rmix(spl, x, v3, dv3, d2r)
            call locate(spl, x, cell, xl)
            if (nq == 1) then
                call evaluate_batch_spline_3d_scalar_quintic_der2_rmix(spl, x, vs, &
                                                                       dvs, d2s)
                call expect3(spl, 1, cell, xl, [0, 0, 0], vs, nbad)
                do j = 1, 3
                    call expect3(spl, 1, cell, xl, unit(j), dvs(j), nbad)
                    call expect3(spl, 1, cell, xl, hessian(:, j), d2s(j), nbad)
                end do
            end if
            do iq = 1, nq
                call expect3(spl, iq, cell, xl, [0, 0, 0], v0(iq), nbad)
                call expect3(spl, iq, cell, xl, [0, 0, 0], v(iq), nbad)
                call expect3(spl, iq, cell, xl, [0, 0, 0], v2(iq), nbad)
                call expect3(spl, iq, cell, xl, [0, 0, 0], v3(iq), nbad)
                do j = 1, 3
                    call expect3(spl, iq, cell, xl, unit(j), dv(j, iq), nbad)
                    call expect3(spl, iq, cell, xl, unit(j), dv2(j, iq), nbad)
                    call expect3(spl, iq, cell, xl, unit(j), dv3(j, iq), nbad)
                    call expect3(spl, iq, cell, xl, hessian(:, j), d2r(j, iq), nbad)
                end do
                do j = 1, 6
                    call expect3(spl, iq, cell, xl, hessian(:, j), d2v(j, iq), nbad)
                end do
            end do
        end do
        if (nbad > 0) then
            write (error_unit, "(a,i0,a,i0,a)") "  3D nq = ", nq, ": ", nbad, &
                " values disagree with the monomial oracle"
            nfail = nfail + 1
        end if
        call destroy_batch_splines_3d(spl)
    end subroutine check_3d

    subroutine check_1d(nq, periodic, nfail)
        integer, intent(in) :: nq
        logical, intent(in) :: periodic
        integer, intent(inout) :: nfail

        integer, parameter :: ngrid = 17, npts = 200
        real(dp), parameter :: x_min = -0.5_dp, x_max = x_min + (ngrid - 1)*0.125_dp
        type(BatchSplineData1D) :: spl
        real(dp) :: y(ngrid, nq), u, x, h, xl, xn
        real(dp) :: v(nq), d1(nq), d2(nq), w(nq), e1(nq), e2(nq), e3(nq)
        real(dp) :: p(2), dp1(2), dp2(2)
        integer :: ip, iq, idx, nbad, j

        call random_seed(put=[(17*j + nq, j=1, 64)])
        call random_number(y)
        if (periodic) y(ngrid, :) = y(1, :)
        call construct_batch_splines_1d(x_min, x_max, y, 5, periodic, spl)
        h = spl%h_step
        nbad = 0
        do ip = 1, npts
            call random_number(u)
            x = x_min + dyadic(u)*(x_max - x_min)
            call evaluate_batch_splines_1d_der2(spl, x, v, d1, d2)
            call evaluate_batch_splines_1d_der3(spl, x, w, e1, e2, e3)
            if (nq == 1) then
                call evaluate_batch_spline_1d_scalar_quintic_der2(spl, x, p(1), &
                                                                  dp1(1), dp2(1))
            else if (nq == 2) then
                call evaluate_batch_spline_1d_pair_quintic_der2(spl, x, p, dp1, dp2)
            end if
            xn = (x - x_min)/h
            idx = max(0, min(ngrid - 2, int(xn)))
            xl = (xn - real(idx, dp))*h
            do iq = 1, nq
                call expect1(spl%coeff(iq, :, idx + 1), xl, 0, v(iq), nbad)
                call expect1(spl%coeff(iq, :, idx + 1), xl, 1, d1(iq), nbad)
                call expect1(spl%coeff(iq, :, idx + 1), xl, 2, d2(iq), nbad)
                call expect1(spl%coeff(iq, :, idx + 1), xl, 0, w(iq), nbad)
                call expect1(spl%coeff(iq, :, idx + 1), xl, 1, e1(iq), nbad)
                call expect1(spl%coeff(iq, :, idx + 1), xl, 2, e2(iq), nbad)
                call expect1(spl%coeff(iq, :, idx + 1), xl, 3, e3(iq), nbad)
                if (nq <= 2) then
                    call expect1(spl%coeff(iq, :, idx + 1), xl, 0, p(iq), nbad)
                    call expect1(spl%coeff(iq, :, idx + 1), xl, 1, dp1(iq), nbad)
                    call expect1(spl%coeff(iq, :, idx + 1), xl, 2, dp2(iq), nbad)
                end if
            end do
        end do
        if (nbad > 0) then
            write (error_unit, "(a,i0,a,i0,a)") "  1D nq = ", nq, ": ", nbad, &
                " values disagree with the monomial oracle"
            nfail = nfail + 1
        end if
        call destroy_batch_splines_1d(spl)
    end subroutine check_1d

    ! Grid steps are 1/8 and points are multiples of 2**-20 of the domain, so
    ! locating a point is exact however the arithmetic is reassociated.
    elemental real(dp) function dyadic(u)
        real(dp), intent(in) :: u
        dyadic = real(int(u*2.0_dp**20), dp)*2.0_dp**(-20)
    end function dyadic

    pure function unit(j) result(a)
        integer, intent(in) :: j
        integer :: a(3)
        a = 0
        a(j) = 1
    end function unit

    ! Interior points only, so no periodic wrap.
    subroutine locate(spl, x, cell, xl)
        type(BatchSplineData3D), intent(in) :: spl
        real(dp), intent(in) :: x(3)
        integer, intent(out) :: cell(3)
        real(dp), intent(out) :: xl(3)
        real(dp) :: xn
        integer :: j

        do j = 1, 3
            xn = (x(j) - spl%x_min(j))*spl%inv_h_step(j)
            cell(j) = max(0, min(spl%num_points(j) - 2, int(xn)))
            xl(j) = (xn - real(cell(j), dp))*spl%h_step(j)
        end do
        cell = cell + 1
    end subroutine locate

    real(qp) function falling(k, a, x) result(f)
        integer, intent(in) :: k, a
        real(dp), intent(in) :: x
        integer :: i

        f = 0
        if (a > k) return
        f = 1
        do i = k - a + 1, k
            f = f*i
        end do
        f = f*real(x, qp)**(k - a)
    end function falling

    subroutine expect3(spl, iq, cell, xl, a, actual, nbad)
        type(BatchSplineData3D), intent(in) :: spl
        integer, intent(in) :: iq, cell(3), a(3)
        real(dp), intent(in) :: xl(3), actual
        integer, intent(inout) :: nbad
        real(qp) :: ref, mag, term
        integer :: k1, k2, k3

        ref = 0
        mag = 0
        do k3 = 0, 5
            do k2 = 0, 5
                do k1 = 0, 5
                    term = real(spl%coeff(iq, k1, k2, k3, cell(1), cell(2), &
                                          cell(3)), qp)*falling(k1, a(1), xl(1))* &
                           falling(k2, a(2), xl(2))*falling(k3, a(3), xl(3))
                    ref = ref + term
                    mag = mag + abs(term)
                end do
            end do
        end do
        call judge(actual, ref, mag, nbad)
    end subroutine expect3

    subroutine expect1(c, xl, a, actual, nbad)
        real(dp), intent(in) :: c(0:5), xl, actual
        integer, intent(in) :: a
        integer, intent(inout) :: nbad
        real(qp) :: ref, mag, term
        integer :: k

        ref = 0
        mag = 0
        do k = 0, 5
            term = real(c(k), qp)*falling(k, a, xl)
            ref = ref + term
            mag = mag + abs(term)
        end do
        call judge(actual, ref, mag, nbad)
    end subroutine expect1

    subroutine judge(actual, ref, mag, nbad)
        real(dp), intent(in) :: actual
        real(qp), intent(in) :: ref, mag
        integer, intent(inout) :: nbad

        if (abs(real(actual, qp) - ref) > 64*epsilon(1.0_dp)*max(mag, &
                                                                real(tiny(1.0_dp), qp))) &
            nbad = nbad + 1
    end subroutine judge
end program test_batch_interpolate_horner
