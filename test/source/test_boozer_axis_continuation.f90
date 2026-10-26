!> Independent power-sum oracle for scalar-primary axis continuation.
program test_boozer_axis_continuation
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use boozer_hybrid, only: boozer_axis_continuation_weights, hybrid_axis_row, &
                             hybrid_radial_row
    implicit none
    integer, parameter :: NS = 17, MODES(4) = [2, 3, 4, 7]
    real(dp), parameter :: XS(5) = [1.0e-200_dp, 1.0e-12_dp, 0.03_dp, &
                                   0.4_dp, 1.0_dp]
    real(dp), parameter :: COEFFICIENTS(3) = [2.0_dp, -3.0_dp, 5.0_dp]
    real(dp) :: w(3, 3), endpoint(3), actual(3), expected(3), identity(3, 3)
    real(dp) :: rv(0:NS - 1), rd(0:NS - 1), rdd(0:NS - 1)
    real(dp) :: v0(0:NS - 1), d0(0:NS - 1), dd0(0:NS - 1)
    real(dp) :: vp(0:NS - 1), dpv(0:NS - 1), ddp(0:NS - 1)
    real(dp) :: vm(0:NS - 1), dmv(0:NS - 1), ddm(0:NS - 1)
    real(dp) :: data(0:NS - 1), direction(0:NS - 1), adjoint(0:NS - 1)
    real(dp) :: cotangent(3), tangent(3), perturbed(3), a, s, x, delta, scale
    integer :: failures, im, ix, ih, j, m, parity
    logical :: half

    failures = 0
    identity = 0
    do j = 1, 3
        identity(j, j) = 1
    end do
    do im = 1, size(MODES)
        m = MODES(im)
        call power_jet(m, 1.0_dp, endpoint)
        call boozer_axis_continuation_weights(m, 1.0_dp, w)
        call check('C2 endpoint weights', maxval(abs(w - identity)), 2.0e-14_dp)
        do ix = 1, size(XS)
            x = XS(ix)
            call boozer_axis_continuation_weights(m, x, w)
            call check('finite weights at small positive s', &
                       merge(0.0_dp, 1.0_dp, all(ieee_is_finite(w))), 0.0_dp)
            actual = matmul(w, endpoint)
            call power_jet(m, x, expected)
            do j = 1, 3
                scale = max(abs(expected(j)), tiny(1.0_dp))
                call check('analytic power-sum jet', &
                           abs(actual(j) - expected(j))/scale, 3.0e-13_dp)
            end do
        end do
        do ih = 0, 1
            half = ih == 1
            a = 1.0_dp/real(NS - 1, dp)
            if (half) a = 0.5_dp*a
            parity = mod(m, 2)
            do j = 0, NS - 1
                s = real(j, dp)/real(NS - 1, dp)
                if (half) s = (real(j, dp) - 0.5_dp)/real(NS - 1, dp)
                data(j) = 0
                if (s > 0) then
                    call power_jet(m, s, expected)
                    data(j) = expected(1)
                end if
                direction(j) = sin(real(j + 1, dp))
            end do
            do ix = 1, 3
                s = a*XS(ix + 2)
                call hybrid_axis_row(NS, s, half, m, rv, rd, rdd)
                actual = [dot_product(rv, data), dot_product(rd, data), &
                          dot_product(rdd, data)]
                call power_jet(m, s, expected)
                do j = 1, 3
                    call check('primary-row analytic power sum', &
                               abs(actual(j)/expected(j) - 1), 2.0e-9_dp)
                end do
            end do
            call hybrid_axis_row(NS, 1.3_dp*a, half, m, rv, rd, rdd)
            call hybrid_radial_row(NS, 1.3_dp*a, half, parity, v0, d0, dd0)
            call check('unchanged rows above join', &
                       maxval(abs(rv - v0)) + maxval(abs(rd - d0)) + &
                       maxval(abs(rdd - dd0)), 0.0_dp)
            call hybrid_axis_row(NS, a, half, m, v0, d0, dd0)
            call hybrid_axis_row(NS, a*(1 - 1.0e-9_dp), half, m, rv, rd, rdd)
            scale = maxval(abs(v0)) + a*maxval(abs(d0)) + &
                    a**2*maxval(abs(dd0))
            call check('C2 row join from below', &
                       (maxval(abs(rv - v0)) + a*maxval(abs(rd - d0)) + &
                        a**2*maxval(abs(rdd - dd0)))/scale, 2.0e-7_dp)

            s = 0.37_dp*a
            delta = 1.0e-5_dp*a
            call hybrid_axis_row(NS, s, half, m, rv, rd, rdd)
            call hybrid_axis_row(NS, s + delta, half, m, vp, dpv, ddp)
            call hybrid_axis_row(NS, s - delta, half, m, vm, dmv, ddm)
            call check('row first derivative by finite difference', &
                       maxval(abs((vp - vm)/(2*delta) - rd))/ &
                       maxval(abs(rd)), 2.0e-8_dp)
            call check('row second derivative by finite difference', &
                       maxval(abs((dpv - dmv)/(2*delta) - rdd))/ &
                       maxval(abs(rdd)), 2.0e-8_dp)
            cotangent = [0.7_dp, -0.3_dp*a, 0.2_dp*a**2]
            adjoint = cotangent(1)*rv + cotangent(2)*rd + cotangent(3)*rdd
            tangent = [dot_product(rv, direction), dot_product(rd, direction), &
                       dot_product(rdd, direction)]
            call check('linear row transpose', &
                       abs(dot_product(cotangent, tangent) - &
                           dot_product(adjoint, direction)), 2.0e-14_dp)
            delta = 1.0e-5_dp
            perturbed = [dot_product(rv, data + delta*direction), &
                         dot_product(rd, data + delta*direction), &
                         dot_product(rdd, data + delta*direction)]
            actual = [dot_product(rv, data - delta*direction), &
                      dot_product(rd, data - delta*direction), &
                      dot_product(rdd, data - delta*direction)]
            call check('primary-data directional finite difference', &
                       abs(dot_product(cotangent, (perturbed - actual)/ &
                           (2*delta) - tangent)), 2.0e-10_dp)
        end do
    end do
    do m = 0, 1
        do ih = 0, 1
            half = ih == 1
            call hybrid_axis_row(NS, 1.0e-5_dp, half, m, rv, rd, rdd)
            call hybrid_radial_row(NS, 1.0e-5_dp, half, m, v0, d0, dd0)
            call check('unchanged m=0,1 rows', &
                       maxval(abs(rv - v0)) + maxval(abs(rd - d0)) + &
                       maxval(abs(rdd - dd0)), 0.0_dp)
        end do
    end do
    if (failures /= 0) error stop 'test_boozer_axis_continuation FAILED'
    print '(a)', 'test_boozer_axis_continuation: all checks passed'

contains

    !> Three independent powers supply the analytic value and s derivatives.
    subroutine power_jet(mode, position, jet)
        integer, intent(in) :: mode
        real(dp), intent(in) :: position
        real(dp), intent(out) :: jet(3)
        real(dp) :: exponent
        integer :: k

        jet = 0
        do k = 1, 3
            exponent = 0.5_dp*real(mode, dp) + k - 1
            jet(1) = jet(1) + COEFFICIENTS(k)*position**exponent
            jet(2) = jet(2) + COEFFICIENTS(k)*exponent*position**(exponent - 1)
            if (mode + 2*(k - 1) /= 2) then
                jet(3) = jet(3) + COEFFICIENTS(k)*exponent* &
                         (exponent - 1)*position**(exponent - 2)
            end if
        end do
    end subroutine power_jet

    subroutine check(label, error, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: error, tolerance

        if (.not. (error <= tolerance)) then
            print '(a,2es14.5)', 'FAIL '//label//': ', error, tolerance
            failures = failures + 1
        end if
    end subroutine check
end program test_boozer_axis_continuation
