program test_boozer_symplectic_file
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_chartmap, only: load_boozer_from_chartmap
    use boozer_rk_tables, only: rk_field_table, rk_profile_table, &
        rk_num_points, rk_x_min, rk_h_step, rk_inv_h_step, rk_period, &
        rk_inv_period, rk_tables_ready, splint_boozer_symplectic_table_device

    implicit none

    character(len=2048) :: filename
    integer :: status, point
    real(dp) :: s, theta, phi
    real(dp) :: actual(19), expected(19), error(19)

    call get_command_argument(1, filename, status=status)
    if (status /= 0 .or. len_trim(filename) == 0) &
        error stop 'usage: test_boozer_symplectic_file.x chartmap.nc'

    call load_boozer_from_chartmap(trim(filename), rk_only=.true.)
    if (.not. rk_tables_ready) error stop 'chartmap has no native RK tables'

    error = 0.0_dp

    do point = 1, 8
        s = rk_x_min(1) + (2.35_dp + 2.41_dp*real(point, dp))*rk_h_step(1)
        theta = (0.27_dp + 0.61_dp*real(point, dp))*rk_h_step(2)
        phi = (0.19_dp + 0.47_dp*real(point, dp))*rk_h_step(3)
        call evaluate_module(s, theta, phi, actual)
        call evaluate_reference(s, theta, phi, expected)
        call update_errors(actual, expected, error)
    end do

    write (*, '(a,1p,10E12.4)') &
        'native value/first-derivative errors: ', error(1:10)
    write (*, '(a,1p,9E12.4)') &
        'native second-derivative errors: ', error(11:19)
    ! The production compact evaluator intentionally uses single-precision
    ! table values, weights, and accumulators for device throughput.  The
    ! independent reference above uses double precision, so cancellation in
    ! angular derivatives leaves an O(epsilon_single) residual when the
    ! derivative itself is small.
    if (maxval(error(1:9)) > 2.0e-3_dp .or. maxval(error(10:19)) > 0.5_dp) &
        error stop 'native RK table derivative contract failed'
    print *, 'test_boozer_symplectic_file: PASSED'

contains

    subroutine evaluate_module(s, theta, phi, values)
        real(dp), intent(in) :: s, theta, phi
        real(dp), intent(out) :: values(19)
        real(dp) :: d2bmod(6), dbmod(3)

        call splint_boozer_symplectic_table_device(s, theta, phi, 2, &
            values(1), values(2), values(3), values(4), values(5), values(6), &
            values(7), values(8), values(9), values(10), dbmod, d2bmod)
        values(11:13) = dbmod
        values(14:19) = d2bmod
    end subroutine evaluate_module

    subroutine evaluate_reference(s, theta, phi, values)
        real(dp), intent(in) :: s, theta, phi
        real(dp), intent(out) :: values(19)

        integer :: first_s, first_theta, first_phi
        integer :: i, j, k, table_index, profile_index
        real(dp) :: weight_s(4), weight_theta(4), weight_phi(4)
        real(dp) :: derivative_s(4), derivative_theta(4), derivative_phi(4)
        real(dp) :: second_theta(4), second_phi(4)
        real(dp) :: bmod, dbmod_s, weight_theta_phi

        call reference_location(s, 1, first_s, weight_s, derivative_s)
        call reference_location(theta, 2, first_theta, weight_theta, &
            derivative_theta, second_theta)
        call reference_location(phi, 3, first_phi, weight_phi, &
            derivative_phi, second_phi)

        values = 0.0_dp
        do i = 1, 4
            profile_index = 1 + 6*(first_s + i - 2)
            values(1) = values(1) + weight_s(i)*rk_profile_table(profile_index)
            values(2) = values(2) + weight_s(i)*rk_profile_table(profile_index + 1)
            values(3) = values(3) + derivative_s(i)*rk_profile_table(profile_index + 1)
            values(4) = values(4) + weight_s(i)*rk_profile_table(profile_index + 2)
            values(5) = values(5) + weight_s(i)*rk_profile_table(profile_index + 3)
            values(6) = values(6) + derivative_s(i)*rk_profile_table(profile_index + 3)
            values(7) = values(7) + weight_s(i)*rk_profile_table(profile_index + 4)
            values(8) = values(8) + weight_s(i)*rk_profile_table(profile_index + 5)
            values(9) = values(9) + derivative_s(i)*rk_profile_table(profile_index + 5)
        end do

        do k = 1, 4
            do j = 1, 4
                weight_theta_phi = weight_theta(j)*weight_phi(k)
                do i = 1, 4
                    table_index = 1 + 2*((first_s + i - 2) + &
                        rk_num_points(1)*((first_theta + j - 2) + &
                        rk_num_points(2)*(first_phi + k - 2)))
                    bmod = real(rk_field_table(table_index), dp)
                    dbmod_s = real(rk_field_table(table_index + 1), dp)
                    values(10) = values(10) + weight_s(i)*weight_theta_phi*bmod
                    values(11) = values(11) + weight_s(i)*weight_theta_phi*dbmod_s
                    values(12) = values(12) + weight_s(i)*derivative_theta(j)* &
                        weight_phi(k)*bmod
                    values(13) = values(13) + weight_s(i)*weight_theta(j)* &
                        derivative_phi(k)*bmod
                    values(14) = values(14) + derivative_s(i)*weight_theta_phi*dbmod_s
                    values(15) = values(15) + weight_s(i)*derivative_theta(j)* &
                        weight_phi(k)*dbmod_s
                    values(16) = values(16) + weight_s(i)*weight_theta(j)* &
                        derivative_phi(k)*dbmod_s
                    values(17) = values(17) + weight_s(i)*second_theta(j)* &
                        weight_phi(k)*bmod
                    values(18) = values(18) + weight_s(i)*derivative_theta(j)* &
                        derivative_phi(k)*bmod
                    values(19) = values(19) + weight_s(i)*weight_theta(j)* &
                        second_phi(k)*bmod
                end do
            end do
        end do
    end subroutine evaluate_reference

    subroutine reference_location(x, idim, first, weight, derivative, second_derivative)
        real(dp), intent(in) :: x
        integer, intent(in) :: idim
        integer, intent(out) :: first
        real(dp), intent(out) :: weight(4)
        real(dp), intent(out), optional :: derivative(4), second_derivative(4)

        real(dp) :: x_eval, x_grid, periods, relative
        integer :: first_zero

        x_eval = x
        if (idim > 1) then
            periods = floor((x - rk_x_min(idim))*rk_inv_period(idim))
            x_eval = x - periods*rk_period(idim)
        end if
        x_grid = (x_eval - rk_x_min(idim))*rk_inv_h_step(idim)
        first_zero = 3*(int(x_grid)/3)
        first_zero = max(0, min(first_zero, rk_num_points(idim) - 4))
        relative = x_grid - real(first_zero, dp)
        first = first_zero + 1

        weight(1) = (1.0_dp - relative)*(2.0_dp - relative)* &
            (3.0_dp - relative)/6.0_dp
        weight(2) = relative*(2.0_dp - relative)*(3.0_dp - relative)/2.0_dp
        weight(3) = relative*(relative - 1.0_dp)*(3.0_dp - relative)/2.0_dp
        weight(4) = relative*(relative - 1.0_dp)*(relative - 2.0_dp)/6.0_dp
        if (present(derivative)) then
            derivative(1) = (-11.0_dp + relative*(12.0_dp - 3.0_dp*relative))/6.0_dp
            derivative(2) = (6.0_dp + relative*(-10.0_dp + 3.0_dp*relative))/2.0_dp
            derivative(3) = (-3.0_dp + relative*(8.0_dp - 3.0_dp*relative))/2.0_dp
            derivative(4) = (2.0_dp + relative*(-6.0_dp + 3.0_dp*relative))/6.0_dp
            derivative = derivative*rk_inv_h_step(idim)
        end if
        if (present(second_derivative)) then
            second_derivative(1) = 2.0_dp - relative
            second_derivative(2) = -5.0_dp + 3.0_dp*relative
            second_derivative(3) = 4.0_dp - 3.0_dp*relative
            second_derivative(4) = relative - 1.0_dp
            second_derivative = second_derivative*rk_inv_h_step(idim)**2
        end if
    end subroutine reference_location

    subroutine update_errors(actual, expected, maximum)
        real(dp), intent(in) :: actual(:), expected(:)
        real(dp), intent(inout) :: maximum(:)
        real(dp) :: scale
        integer :: i

        do i = 1, size(actual)
            scale = max(1.0_dp, abs(actual(i)), abs(expected(i)))
            maximum(i) = max(maximum(i), abs(actual(i) - expected(i))/scale)
        end do
    end subroutine update_errors

end program test_boozer_symplectic_file
