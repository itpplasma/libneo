program test_boozer_symplectic_table
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
    use boozer_rk_tables, only: rk_field_table, rk_profile_table, &
        rk_num_points, rk_x_min, rk_h_step, rk_inv_h_step, rk_period, &
        rk_inv_period, rk_tables_ready, splint_boozer_symplectic_table_device

    implicit none

    integer, parameter :: n = 7
    real(dp), parameter :: eps = 2.0e-3_dp
    real(dp) :: base(19), plus(19), minus(19), finite_difference
    integer :: i, j, k, point, offset

    allocate (rk_field_table(2*n*n*n), rk_profile_table(6*n))
    do k = 1, n
        do j = 1, n
            do i = 1, n
                point = (i - 1) + n*((j - 1) + n*(k - 1))
                offset = 1 + 2*point
                rk_field_table(offset) = real(2.0_dp + 0.03_dp*i*j + &
                    0.02_dp*j*k + 0.01_dp*i*k, sp)
                rk_field_table(offset + 1) = real(-0.4_dp + 0.02_dp*i*j - &
                    0.015_dp*j*k + 0.025_dp*i*k, sp)
            end do
        end do
    end do
    do i = 1, n
        offset = 1 + 6*(i - 1)
        rk_profile_table(offset:offset + 5) = real([ &
            0.2_dp*i*i, 0.1_dp*i*i + 0.03_dp*i, &
            -0.4_dp*i, 0.07_dp*i*i - 0.02_dp*i, &
            1.5_dp + 0.2_dp*i, -0.05_dp*i*i + 0.04_dp*i], sp)
    end do
    rk_num_points = n
    rk_x_min = [0.0_dp, 0.0_dp, 0.0_dp]
    rk_h_step = [0.2_dp, 0.3_dp, 0.25_dp]
    rk_inv_h_step = 1.0_dp/rk_h_step
    rk_period = rk_h_step*real(rk_num_points - 1, dp)
    rk_inv_period = 1.0_dp/rk_period
    rk_tables_ready = .true.

    call evaluate(0.62_dp, 0.43_dp, 0.77_dp, 2, base)

    call evaluate(0.62_dp + eps, 0.43_dp, 0.77_dp, 2, plus)
    call evaluate(0.62_dp - eps, 0.43_dp, 0.77_dp, 2, minus)
    finite_difference = (plus(5) - minus(5))/(2.0_dp*eps)
    call require_close('d2B/ds2', base(8), finite_difference)
    call require_close('d2Aphi/ds2', base(3), &
        (plus(2) - minus(2))/(2.0_dp*eps))
    call require_close('d2Btheta/ds2', base(16), &
        (plus(15) - minus(15))/(2.0_dp*eps))
    call require_close('d2Bphi/ds2', base(19), &
        (plus(18) - minus(18))/(2.0_dp*eps))

    call evaluate(0.62_dp, 0.43_dp + eps, 0.77_dp, 2, plus)
    call evaluate(0.62_dp, 0.43_dp - eps, 0.77_dp, 2, minus)
    call require_close('d2B/dsdtheta', base(9), &
        (plus(5) - minus(5))/(2.0_dp*eps))
    call require_close('d2B/dtheta2', base(11), &
        (plus(6) - minus(6))/(2.0_dp*eps))
    call require_close('d2B/dthetadphi theta', base(12), &
        (plus(7) - minus(7))/(2.0_dp*eps))

    call evaluate(0.62_dp, 0.43_dp, 0.77_dp + eps, 2, plus)
    call evaluate(0.62_dp, 0.43_dp, 0.77_dp - eps, 2, minus)
    call require_close('d2B/dsdphi', base(10), &
        (plus(5) - minus(5))/(2.0_dp*eps))
    call require_close('d2B/dthetadphi phi', base(12), &
        (plus(6) - minus(6))/(2.0_dp*eps))
    call require_close('d2B/dphi2', base(13), &
        (plus(7) - minus(7))/(2.0_dp*eps))

    call evaluate(0.62_dp, 0.43_dp, 0.77_dp, 3, plus)
    if (any(abs(plus(11:13)) > epsilon(1.0_sp))) &
        error stop 'radial-mixed mode computed angular Hessian entries'
    call evaluate(0.62_dp, 0.43_dp, 0.77_dp, 0, plus)
    if (abs(plus(3)) > epsilon(1.0_sp) .or. &
            any(abs(plus(8:13)) > epsilon(1.0_sp)) .or. &
            abs(plus(16)) > epsilon(1.0_sp) .or. &
            abs(plus(19)) > epsilon(1.0_sp)) &
        error stop 'first-derivative mode computed Hessian entries'

    print *, 'compact symplectic Boozer derivative oracle passed'

contains

    subroutine evaluate(s, theta, phi, mode_secders, values)
        real(dp), intent(in) :: s, theta, phi
        integer, intent(in) :: mode_secders
        real(dp), intent(out) :: values(19)
        real(dp) :: d2bmod(6)

        call splint_boozer_symplectic_table_device(s, theta, phi, &
            mode_secders, values(1), values(2), values(3), values(14), &
            values(15), values(16), values(17), values(18), values(19), &
            values(4), values(5:7), d2bmod)
        values(8:13) = d2bmod
    end subroutine evaluate

    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: scale

        scale = max(1.0_dp, abs(actual), abs(expected))
        if (abs(actual - expected) > 2.0e-3_dp*scale) then
            print *, trim(label), ': ', actual, expected
            error stop 'compact symplectic table derivative mismatch'
        end if
    end subroutine require_close

end program test_boozer_symplectic_table
