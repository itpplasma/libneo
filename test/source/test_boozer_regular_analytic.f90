!> Analytic oracles for the regular near-axis Boozer field (boozer_regular).
!>
!> The grid data are exact rho**|m| harmonics with polynomial profiles in s.
!> Checks: value, gradient and Hessian at X = Y = 0 against the analytic
!> Taylor coefficients (1e-10); spline accuracy off the axis; derivatives
!> against finite differences; theta-independent limits through the axis; and
!> field-period periodicity for nfp = 1, 2, 5 and both helicity signs.
program test_boozer_regular_analytic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use interpolate, only: BatchSplineData1D, construct_batch_splines_1d
    use boozer_regular, only: boozer_regular_t, regular_field_t, &
                              build_boozer_regular, REGULAR_VALUE, &
                              REGULAR_FIRST, REGULAR_SECOND
    implicit none

    real(dp), parameter :: TWOPI = 6.283185307179586476925286766559_dp
    real(dp), parameter :: TORFLUX = 0.5_dp
    real(dp), parameter :: ALPHA = 0.7_dp, BETA = -1.1_dp
    integer, parameter :: N_RHO = 41, N_THETA = 33, N_PHI = 25
    integer :: n_fail

    n_fail = 0
    call test_axis_taylor()
    call test_off_axis_accuracy()
    call test_derivatives_finite_difference()
    call test_axis_limits()
    call test_field_period(1, 1)
    call test_field_period(2, 1)
    call test_field_period(5, 1)
    call test_field_period(1, -1)
    call test_field_period(2, -1)
    call test_field_period(5, -1)
    if (n_fail > 0) then
        print *, 'test_boozer_regular_analytic: FAILED checks:', n_fail
        error stop 1
    end if
    print *, 'test_boozer_regular_analytic: all checks passed'

contains

    !> |B| with exact m = 0..3 near-axis structure; helical terms when eps /= 0.
    pure function bmod_exact(rho, theta, phi, nfp, sgn, eps) result(b)
        real(dp), intent(in) :: rho, theta, phi, eps
        integer, intent(in) :: nfp, sgn
        real(dp) :: b, s, zeta

        s = rho**2
        zeta = real(sgn*nfp, dp)*phi
        b = 1.0_dp + 0.2_dp*s - 0.05_dp*s**2 &
            + 0.3_dp*rho*(1.0_dp + 0.5_dp*s + 0.2_dp*s**2)*cos(theta - ALPHA) &
            + 0.15_dp*s*(1.0_dp - 0.3_dp*s + 0.1_dp*s**2)*cos(2.0_dp*theta - BETA) &
            + 0.05_dp*rho**3*cos(3.0_dp*theta) &
            + eps*(0.1_dp*cos(zeta) + 0.04_dp*rho*cos(theta - zeta) &
                   + 0.03_dp*s*cos(2.0_dp*theta + zeta))
    end function bmod_exact

    pure real(dp) function itor(s)
        real(dp), intent(in) :: s
        itor = s*(0.1_dp + 0.05_dp*s)
    end function itor

    pure real(dp) function gpol(s)
        real(dp), intent(in) :: s
        gpol = 2.0_dp + 0.3_dp*s + 0.1_dp*s**2
    end function gpol

    pure real(dp) function aphi_exact(s)
        real(dp), intent(in) :: s
        aphi_exact = -TORFLUX*(0.4_dp*s + 0.1_dp*s**2)
    end function aphi_exact

    subroutine build_field(field, nfp, sgn, eps)
        type(boozer_regular_t), intent(inout) :: field
        integer, intent(in) :: nfp, sgn
        real(dp), intent(in) :: eps
        real(dp) :: rho(N_RHO), bmod(N_RHO, N_THETA, N_PHI)
        real(dp) :: btheta(N_RHO), bphi(N_RHO), ap(N_RHO, 1)
        type(BatchSplineData1D) :: aphi
        integer :: i, j, k

        do i = 1, N_RHO
            rho(i) = real(i - 1, dp)/real(N_RHO - 1, dp)
            btheta(i) = itor(rho(i)**2)
            bphi(i) = gpol(rho(i)**2)
            ap(i, 1) = aphi_exact(rho(i))
            do k = 1, N_PHI
                do j = 1, N_THETA
                    bmod(i, j, k) = bmod_exact(rho(i), &
                        TWOPI*real(j - 1, dp)/real(N_THETA - 1, dp), &
                        TWOPI*real(k - 1, dp)/real((N_PHI - 1)*nfp, dp), &
                        nfp, sgn, eps)
                end do
            end do
        end do
        ! A_phi is splined over a uniform s grid; rho(i) doubles as s(i).
        call construct_batch_splines_1d(0.0_dp, 1.0_dp, ap, 5, .false., aphi)
        call build_boozer_regular(field, rho, bmod, btheta, bphi, aphi, TORFLUX, &
                                  nfp, 5, 5)
    end subroutine build_field

    subroutine check(label, got, want, tol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, want, tol

        if (.not. (abs(got - want) <= tol*max(1.0_dp, abs(want)) + 1.0e-13_dp)) then
            print '(a,a,3es16.8)', 'FAIL ', label, got, want, abs(got - want)
            n_fail = n_fail + 1
        end if
    end subroutine check

    subroutine check_array(label, got, want, tol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got(:), want(:), tol
        integer :: k
        character(len=8) :: idx

        do k = 1, size(got)
            write (idx, '(i0)') k
            call check(label//'('//trim(idx)//')', got(k), want(k), tol)
        end do
    end subroutine check_array

    !> Quotient q = n/b for a phi-independent numerator; the independent
    !> oracle for h = B_k/|B| (first and second derivatives).
    subroutine quotient(n, nd, nd2, b, bd, bd2, q, qd, qd2)
        real(dp), intent(in) :: n, nd(3), nd2(6), b, bd(3), bd2(6)
        real(dp), intent(out) :: q, qd(3), qd2(6)
        integer, parameter :: pi(6) = [1, 1, 1, 2, 2, 3], pj(6) = [1, 2, 3, 2, 3, 3]
        integer :: k

        q = n/b
        qd = nd/b - n*bd/b**2
        do k = 1, 6
            qd2(k) = nd2(k)/b - (nd(pi(k))*bd(pj(k)) + nd(pj(k))*bd(pi(k)))/b**2 &
                     - n*bd2(k)/b**2 + 2.0_dp*n*bd(pi(k))*bd(pj(k))/b**3
        end do
    end subroutine quotient

    subroutine test_axis_taylor()
        type(boozer_regular_t) :: field
        type(regular_field_t) :: f
        real(dp) :: b, bd(3), bd2(6), q, qd(3), qd2(6)
        real(dp), parameter :: tol = 1.0e-10_dp
        integer :: k
        real(dp) :: phi, err

        call build_field(field, 1, 1, 0.0_dp)
        b = 1.0_dp
        bd = [0.3_dp*cos(ALPHA), 0.3_dp*sin(ALPHA), 0.0_dp]
        bd2 = [0.4_dp + 0.3_dp*cos(BETA), 0.3_dp*sin(BETA), 0.0_dp, &
               0.4_dp - 0.3_dp*cos(BETA), 0.0_dp, 0.0_dp]
        err = 0.0_dp
        do k = 0, 3
            phi = 0.37_dp + 1.3_dp*real(k, dp)
            call field%eval(0.0_dp, 0.0_dp, phi, REGULAR_SECOND, f)
            err = max(err, abs(f%Bmod - b), maxval(abs(f%dBmod - bd)), &
                      maxval(abs(f%d2Bmod - bd2)))
            call check('axis B', f%Bmod, b, tol)
            call check_array('axis dB', f%dBmod, bd, tol)
            call check_array('axis d2B', f%d2Bmod, bd2, tol)
            call quotient(0.0_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
                          [0.2_dp, 0.0_dp, 0.0_dp, 0.2_dp, 0.0_dp, 0.0_dp], &
                          b, bd, bd2, q, qd, qd2)
            call check('axis hth', f%hth, q, tol)
            call check_array('axis dhth', f%dhth, qd, tol)
            call check_array('axis d2hth', f%d2hth, qd2, tol)
            call quotient(2.0_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
                          [0.6_dp, 0.0_dp, 0.0_dp, 0.6_dp, 0.0_dp, 0.0_dp], &
                          b, bd, bd2, q, qd, qd2)
            call check('axis hph', f%hph, q, tol)
            call check_array('axis dhph', f%dhph, qd, tol)
            call check_array('axis d2hph', f%d2hph, qd2, tol)
            call check('axis Ath', f%Ath, 0.0_dp, tol)
            call check_array('axis dAth', f%dAth, [0.0_dp, 0.0_dp, 0.0_dp], tol)
            call check_array('axis d2Ath', f%d2Ath, &
                             [2*TORFLUX, 0.0_dp, 0.0_dp, 2*TORFLUX, 0.0_dp, 0.0_dp], &
                             tol)
            call check('axis Aph', f%Aph, 0.0_dp, tol)
            call check_array('axis dAph', f%dAph, [0.0_dp, 0.0_dp, 0.0_dp], tol)
            call check_array('axis d2Aph', f%d2Aph, &
                             [-0.8_dp*TORFLUX, 0.0_dp, 0.0_dp, -0.8_dp*TORFLUX, &
                              0.0_dp, 0.0_dp], tol)
        end do
        print '(a,es10.2)', 'axis max error of |B| and its derivatives:', err
    end subroutine test_axis_taylor

    subroutine test_off_axis_accuracy()
        type(boozer_regular_t) :: field
        type(regular_field_t) :: f
        real(dp) :: rho, theta, phi, err_b, err_h
        integer :: k

        call build_field(field, 2, 1, 1.0_dp)
        err_b = 0.0_dp
        err_h = 0.0_dp
        do k = 1, 200
            rho = 0.9_dp*real(k, dp)/200.0_dp
            theta = 0.913_dp*real(k, dp)
            phi = 0.471_dp*real(k, dp)
            call field%eval(rho*cos(theta), rho*sin(theta), phi, REGULAR_VALUE, f)
            err_b = max(err_b, abs(f%Bmod - bmod_exact(rho, theta, phi, 2, 1, &
                                                       1.0_dp)))
            err_h = max(err_h, abs(f%hph - gpol(rho**2)/bmod_exact(rho, theta, &
                                                                   phi, 2, 1, 1.0_dp)))
            err_h = max(err_h, abs(f%hth - itor(rho**2)/bmod_exact(rho, theta, &
                                                                   phi, 2, 1, 1.0_dp)))
        end do
        print '(a,2es10.2)', 'off-axis max error |B|, h:', err_b, err_h
        call check('off-axis |B| spline error', err_b, 0.0_dp, 1.0e-6_dp)
        call check('off-axis h spline error', err_h, 0.0_dp, 1.0e-6_dp)
    end subroutine test_off_axis_accuracy

    !> First derivatives against central differences of the values, second
    !> derivatives against central differences of the first derivatives.
    subroutine test_derivatives_finite_difference()
        type(boozer_regular_t) :: field
        type(regular_field_t) :: f, fp, fm
        real(dp), parameter :: dh = 1.0e-6_dp
        real(dp) :: p(3), e(3), d1(5), d2(5, 3)
        integer :: ip, dir, k
        integer, parameter :: col(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
        real(dp), parameter :: points(3, 4) = reshape([ &
            0.3_dp, -0.2_dp, 0.4_dp, 1.0e-3_dp, 2.0e-3_dp, 1.0_dp, &
            -0.5_dp, 0.1_dp, 2.5_dp, 0.0_dp, 0.0_dp, 0.2_dp], [3, 4])

        call build_field(field, 5, -1, 1.0_dp)
        do ip = 1, 4
            p = points(:, ip)
            call field%eval(p(1), p(2), p(3), REGULAR_SECOND, f)
            do dir = 1, 3
                e = 0.0_dp
                e(dir) = dh
                call field%eval(p(1) + e(1), p(2) + e(2), p(3) + e(3), &
                                REGULAR_FIRST, fp)
                call field%eval(p(1) - e(1), p(2) - e(2), p(3) - e(3), &
                                REGULAR_FIRST, fm)
                d1 = [fp%Bmod - fm%Bmod, fp%hth - fm%hth, fp%hph - fm%hph, &
                      fp%Ath - fm%Ath, fp%Aph - fm%Aph]/(2*dh)
                call check('FD dB', f%dBmod(dir), d1(1), 1.0e-7_dp)
                call check('FD dhth', f%dhth(dir), d1(2), 1.0e-7_dp)
                call check('FD dhph', f%dhph(dir), d1(3), 1.0e-7_dp)
                call check('FD dAth', f%dAth(dir), d1(4), 1.0e-7_dp)
                call check('FD dAph', f%dAph(dir), d1(5), 1.0e-7_dp)
                do k = 1, 3
                    d2(:, k) = [fp%dBmod(k) - fm%dBmod(k), fp%dhth(k) - fm%dhth(k), &
                                fp%dhph(k) - fm%dhph(k), fp%dAth(k) - fm%dAth(k), &
                                fp%dAph(k) - fm%dAph(k)]/(2*dh)
                    call check('FD d2B', f%d2Bmod(col(dir, k)), d2(1, k), 1.0e-6_dp)
                    call check('FD d2hth', f%d2hth(col(dir, k)), d2(2, k), 1.0e-6_dp)
                    call check('FD d2hph', f%d2hph(col(dir, k)), d2(3, k), 1.0e-6_dp)
                    call check('FD d2Ath', f%d2Ath(col(dir, k)), d2(4, k), 1.0e-6_dp)
                    call check('FD d2Aph', f%d2Aph(col(dir, k)), d2(5, k), 1.0e-6_dp)
                end do
            end do
        end do
    end subroutine test_derivatives_finite_difference

    !> Values, gradients and Hessians approach the axis values from every
    !> direction: theta-independent limits, gradient continuous to O(rho),
    !> Hessian continuous (rho**3 R is C^2).
    subroutine test_axis_limits()
        type(boozer_regular_t) :: field
        type(regular_field_t) :: f0, f
        real(dp) :: rho, theta
        integer :: k, l

        call build_field(field, 2, 1, 1.0_dp)
        call field%eval(0.0_dp, 0.0_dp, 0.8_dp, REGULAR_SECOND, f0)
        do l = 1, 3
            rho = 10.0_dp**(-3*l)
            do k = 0, 15
                theta = TWOPI*real(k, dp)/16.0_dp + 0.1_dp
                call field%eval(rho*cos(theta), rho*sin(theta), 0.8_dp, &
                                REGULAR_SECOND, f)
                call check('limit B', f%Bmod, f0%Bmod + rho*(cos(theta)*f0%dBmod(1) &
                           + sin(theta)*f0%dBmod(2)), 10*rho**2)
                call check('limit hth', f%hth, f0%hth, 10*rho**2)
                call check('limit hph', f%hph, f0%hph, 10*rho)
                call check('limit Ath', f%Ath, f0%Ath, 10*rho**2)
                call check('limit Aph', f%Aph, f0%Aph, 10*rho**2)
                call check_array('limit dB', f%dBmod, f0%dBmod, 10*rho)
                call check_array('limit dhph', f%dhph, f0%dhph, 10*rho)
                call check_array('limit d2B', f%d2Bmod, f0%d2Bmod, 10*rho)
                call check_array('limit d2hth', f%d2hth, f0%d2hth, 10*rho)
            end do
        end do
    end subroutine test_axis_limits

    !> phi and phi + k 2 pi/nfp give identical values and derivatives, and
    !> the field matches the analytic helical data (wrong period or sign
    !> would fail the latter).
    subroutine test_field_period(nfp, sgn)
        integer, intent(in) :: nfp, sgn
        type(boozer_regular_t) :: field
        type(regular_field_t) :: f, g
        real(dp) :: rho, theta, phi, period, err
        integer :: k, shift
        integer, parameter :: shifts(4) = [-3, -1, 1, 7]

        call build_field(field, nfp, sgn, 1.0_dp)
        period = TWOPI/real(nfp, dp)
        err = 0.0_dp
        do k = 1, 40
            rho = 0.02_dp*real(k, dp)
            theta = 1.7_dp*real(k, dp)
            phi = -2.0_dp + 0.31_dp*real(k, dp)
            if (k == 1) rho = 0.0_dp
            call field%eval(rho*cos(theta), rho*sin(theta), phi, REGULAR_SECOND, f)
            err = max(err, abs(f%Bmod - bmod_exact(rho, theta, phi, nfp, sgn, &
                                                   1.0_dp)))
            do shift = 1, size(shifts)
                call field%eval(rho*cos(theta), rho*sin(theta), &
                                phi + real(shifts(shift), dp)*period, &
                                REGULAR_SECOND, g)
                call check('period B', g%Bmod, f%Bmod, 1.0e-12_dp)
                call check_array('period dB', g%dBmod, f%dBmod, 1.0e-11_dp)
                call check_array('period d2B', g%d2Bmod, f%d2Bmod, 1.0e-10_dp)
                call check_array('period dhph', g%dhph, f%dhph, 1.0e-11_dp)
            end do
        end do
        print '(a,i2,a,i3,a,es10.2)', 'nfp=', nfp, ' sign=', sgn, &
            ' max |B| error vs analytic:', err
        call check('period analytic', err, 0.0_dp, 1.0e-6_dp)
    end subroutine test_field_period

end program test_boozer_regular_analytic
