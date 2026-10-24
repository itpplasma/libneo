!> Independent polynomial oracles for radial second derivatives and their
!> transpose, including the factored odd-parity interpolation near the axis.
program test_sfl_radial_second
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use boozer_hybrid, only: hybrid_radial_row, hybrid_radial_second, &
                             hybrid_radial_second_t
    implicit none

    integer, parameter :: NS = 129
    real(dp), parameter :: HS = 1.0_dp/real(NS - 1, dp)
    real(dp), parameter :: SURFACES(4) = [HS/2, 0.245_dp, 0.495_dp, 1 - HS/2]
    real(dp) :: nodes(0:NS - 1), data(0:NS - 1), tangent(0:NS - 1)
    real(dp) :: rv(0:NS - 1), rd(0:NS - 1), rdd(0:NS - 1), cb(0:NS - 1)
    real(dp) :: rough(0:NS - 1), g(0:NS - 1), surface, value, expected
    real(dp) :: lhs, rhs, scale, error, base
    integer :: j, parity, is, failures

    failures = 0
    nodes = [(real(j, dp)*HS, j=0, NS - 1)]
    rough = [(sin(real(j, dp) + 0.3_dp), j=0, NS - 1)]
    do parity = 0, 1
        g = 1.0_dp
        if (parity == 1) g = sqrt(nodes)
        base = merge(8.0_dp, 0.0625_dp, parity == 0)
        data = g*(base + 0.125_dp*nodes + 0.0625_dp*nodes**2)
        tangent = g*(0.125_dp - 0.5_dp*nodes + 0.25_dp*nodes**2)
        do is = 1, size(SURFACES)
            surface = SURFACES(is)
            call hybrid_radial_row(NS, surface, .false., parity, rv, rd, rdd)
            value = hybrid_radial_second(rdd, data, parity)
            expected = second_polynomial(surface, parity, base, 0.125_dp, 0.0625_dp)
            call check('polynomial second derivative', value, expected, 1.0e-9_dp)
            call check('linear tangent second derivative', &
                       hybrid_radial_second(rdd, tangent, parity), &
                       second_polynomial(surface, parity, 0.125_dp, -0.5_dp, 0.25_dp), &
                       1.0e-9_dp)
            if (parity == 0) then
                call check('constant has zero second derivative', &
                           hybrid_radial_second(rdd, 8.0_dp*g, parity), 0.0_dp, &
                           1.0e-13_dp)
            end if
            cb = 0.0_dp
            call hybrid_radial_second_t(rdd, parity, 0.75_dp, cb)
            lhs = 0.75_dp*hybrid_radial_second(rdd, rough, parity)
            rhs = dot_product(cb, rough)
            scale = max(tiny(1.0_dp), norm2(cb)*norm2(rough))
            error = abs(lhs - rhs)/scale
            call check('transpose duality', error, 0.0_dp, 1.0e-13_dp)
        end do
    end do
    if (failures /= 0) error stop 'test_sfl_radial_second FAILED'
    print '(a)', 'test_sfl_radial_second PASSED'

contains

    pure real(dp) function second_polynomial(s, p, a, b, c) result(value)
        real(dp), intent(in) :: s, a, b, c
        integer, intent(in) :: p

        value = 2*c
        if (p == 1) value = -0.25_dp*a/s**1.5_dp + &
                            0.75_dp*b/sqrt(s) + 3.75_dp*c*sqrt(s)
    end function second_polynomial

    subroutine check(label, actual, reference, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, reference, tolerance
        real(dp) :: relative_error

        relative_error = abs(actual - reference)/max(1.0_dp, abs(reference))
        if (.not. ieee_is_finite(actual) .or. relative_error > tolerance) then
            print '(a, 2i4, f10.6, 3es15.6)', 'FAIL '//label, parity, is, surface, &
                actual, reference, relative_error
            failures = failures + 1
        end if
    end subroutine check

end program test_sfl_radial_second
