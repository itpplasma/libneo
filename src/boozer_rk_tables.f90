module boozer_rk_tables
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32

    implicit none
    private

    public :: rk_field_table, rk_profile_table
    public :: rk_num_points, rk_x_min, rk_h_step, rk_inv_h_step, rk_period
    public :: rk_inv_period
    public :: rk_tables_ready
    public :: splint_boozer_rk_table_device
    public :: splint_boozer_symplectic_table_device

    ! Keep the compact Landreman-compatible tables in a small module of their
    ! own. NVHPC 26.5 otherwise combines these declarations with the generic
    ! spline descriptors and emits an invalid device-side static initializer.
    real(sp), allocatable, save :: rk_field_table(:)
    real(sp), allocatable, save :: rk_profile_table(:)
    integer, save :: rk_num_points(3)
    real(dp), save :: rk_x_min(3)
    real(dp), save :: rk_h_step(3)
    real(dp), save :: rk_inv_h_step(3)
    real(dp), save :: rk_period(3)
    real(dp), save :: rk_inv_period(3)
    logical, save :: rk_tables_ready = .false.
    !$acc declare create(rk_field_table, rk_profile_table, rk_num_points, &
    !$acc& rk_x_min, rk_h_step, rk_inv_h_step, rk_period, rk_inv_period, &
    !$acc& rk_tables_ready)

contains

    subroutine cubic_table_location(x, idim, first, weight, derivative, &
            second_derivative)
        !$acc routine seq
        real(dp), intent(in) :: x
        integer, intent(in) :: idim
        integer, intent(out) :: first
        real(sp), intent(out) :: weight(4)
        real(sp), intent(out), optional :: derivative(4)
        real(sp), intent(out), optional :: second_derivative(4)

        real(dp) :: x_eval, x_grid, periods, relative
        real(sp) :: relative_sp
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
        relative_sp = real(relative, sp)
        first = first_zero + 1

        weight(1) = (1.0_sp - relative_sp)*(2.0_sp - relative_sp)* &
            (3.0_sp - relative_sp)/6.0_sp
        weight(2) = relative_sp*(2.0_sp - relative_sp)* &
            (3.0_sp - relative_sp)/2.0_sp
        weight(3) = relative_sp*(relative_sp - 1.0_sp)* &
            (3.0_sp - relative_sp)/2.0_sp
        weight(4) = relative_sp*(relative_sp - 1.0_sp)* &
            (relative_sp - 2.0_sp)/6.0_sp
        if (present(derivative)) then
            derivative(1) = (-11.0_sp + relative_sp* &
                (12.0_sp - 3.0_sp*relative_sp))/6.0_sp
            derivative(2) = (6.0_sp + relative_sp* &
                (-10.0_sp + 3.0_sp*relative_sp))/2.0_sp
            derivative(3) = (-3.0_sp + relative_sp* &
                (8.0_sp - 3.0_sp*relative_sp))/2.0_sp
            derivative(4) = (2.0_sp + relative_sp* &
                (-6.0_sp + 3.0_sp*relative_sp))/6.0_sp
            derivative = derivative*real(rk_inv_h_step(idim), sp)
        end if
        if (present(second_derivative)) then
            second_derivative(1) = 2.0_sp - relative_sp
            second_derivative(2) = -5.0_sp + 3.0_sp*relative_sp
            second_derivative(3) = 4.0_sp - 3.0_sp*relative_sp
            second_derivative(4) = relative_sp - 1.0_sp
            second_derivative = second_derivative* &
                real(rk_inv_h_step(idim)**2, sp)
        end if
    end subroutine cubic_table_location

    subroutine evaluate_rk_profile_table(s, value)
        !$acc routine seq
        real(dp), intent(in) :: s
        real(sp), intent(out) :: value(6)

        integer :: first, i, iq, table_index
        real(sp) :: weight(4)

        call cubic_table_location(s, 1, first, weight)
        value = 0.0_sp
        do i = 1, 4
            do iq = 1, 6
                table_index = iq + 6*(first + i - 2)
                value(iq) = value(iq) + weight(i)*rk_profile_table(table_index)
            end do
        end do
    end subroutine evaluate_rk_profile_table

    subroutine evaluate_rk_field_table(s, theta, phi, value)
        !$acc routine seq
        real(dp), intent(in) :: s, theta, phi
        real(sp), intent(out) :: value(4)

        integer :: first_s, first_theta, first_phi
        integer :: i, j, k, table_index
        real(sp) :: weight_s(4), weight_theta(4), weight_phi(4)
        real(sp) :: derivative_theta(4), derivative_phi(4)
        real(sp) :: bmod, weight_theta_phi

        call cubic_table_location(s, 1, first_s, weight_s)
        call cubic_table_location(theta, 2, first_theta, weight_theta, &
            derivative_theta)
        call cubic_table_location(phi, 3, first_phi, weight_phi, derivative_phi)
        value = 0.0_sp
        do k = 1, 4
            do j = 1, 4
                weight_theta_phi = weight_theta(j)*weight_phi(k)
                do i = 1, 4
                    table_index = 1 + 2*((first_s + i - 2) + &
                        rk_num_points(1)*((first_theta + j - 2) + &
                        rk_num_points(2)*(first_phi + k - 2)))
                    bmod = rk_field_table(table_index)
                    value(1) = value(1) + &
                        weight_s(i)*weight_theta_phi*bmod
                    value(2) = value(2) + &
                        weight_s(i)*weight_theta_phi* &
                        rk_field_table(table_index + 1)
                    value(3) = value(3) + weight_s(i)* &
                        derivative_theta(j)*weight_phi(k)*bmod
                    value(4) = value(4) + weight_s(i)*weight_theta(j)* &
                        derivative_phi(k)*bmod
                end do
            end do
        end do
    end subroutine evaluate_rk_field_table

    subroutine evaluate_symplectic_profile_table(s, mode_secders, value, &
            second_value)
        !$acc routine seq
        real(dp), intent(in) :: s
        integer, intent(in) :: mode_secders
        real(sp), intent(out) :: value(6), second_value(3)

        integer :: first, i, iq, table_index
        real(sp) :: weight(4), derivative(4)

        if (mode_secders > 0) then
            call cubic_table_location(s, 1, first, weight, derivative)
        else
            call cubic_table_location(s, 1, first, weight)
        end if
        value = 0.0_sp
        second_value = 0.0_sp
        do i = 1, 4
            do iq = 1, 6
                table_index = iq + 6*(first + i - 2)
                value(iq) = value(iq) + weight(i)*rk_profile_table(table_index)
            end do
            if (mode_secders > 0) then
                second_value(1) = second_value(1) + derivative(i)* &
                    rk_profile_table(2 + 6*(first + i - 2))
                second_value(2) = second_value(2) + derivative(i)* &
                    rk_profile_table(4 + 6*(first + i - 2))
                second_value(3) = second_value(3) + derivative(i)* &
                    rk_profile_table(6 + 6*(first + i - 2))
            end if
        end do
    end subroutine evaluate_symplectic_profile_table

    subroutine evaluate_symplectic_field_table(s, theta, phi, mode_secders, &
            value, second_value)
        !$acc routine seq
        real(dp), intent(in) :: s, theta, phi
        integer, intent(in) :: mode_secders
        real(sp), intent(out) :: value(4), second_value(6)

        integer :: first_s, first_theta, first_phi
        integer :: i, j, k, table_index
        real(sp) :: weight_s(4), weight_theta(4), weight_phi(4)
        real(sp) :: derivative_s(4), derivative_theta(4), derivative_phi(4)
        real(sp) :: second_theta(4), second_phi(4)
        real(sp) :: bmod, dbmod_s, weight_theta_phi

        if (mode_secders > 0) then
            call cubic_table_location(s, 1, first_s, weight_s, derivative_s)
        else
            call cubic_table_location(s, 1, first_s, weight_s)
        end if
        if (mode_secders == 2) then
            call cubic_table_location(theta, 2, first_theta, weight_theta, &
                derivative_theta, second_theta)
            call cubic_table_location(phi, 3, first_phi, weight_phi, &
                derivative_phi, second_phi)
        else
            call cubic_table_location(theta, 2, first_theta, weight_theta, &
                derivative_theta)
            call cubic_table_location(phi, 3, first_phi, weight_phi, &
                derivative_phi)
        end if
        value = 0.0_sp
        second_value = 0.0_sp
        do k = 1, 4
            do j = 1, 4
                weight_theta_phi = weight_theta(j)*weight_phi(k)
                do i = 1, 4
                    table_index = 1 + 2*((first_s + i - 2) + &
                        rk_num_points(1)*((first_theta + j - 2) + &
                        rk_num_points(2)*(first_phi + k - 2)))
                    bmod = rk_field_table(table_index)
                    dbmod_s = rk_field_table(table_index + 1)
                    value(1) = value(1) + &
                        weight_s(i)*weight_theta_phi*bmod
                    value(2) = value(2) + &
                        weight_s(i)*weight_theta_phi*dbmod_s
                    value(3) = value(3) + weight_s(i)* &
                        derivative_theta(j)*weight_phi(k)*bmod
                    value(4) = value(4) + weight_s(i)*weight_theta(j)* &
                        derivative_phi(k)*bmod
                    if (mode_secders > 0) then
                        second_value(1) = second_value(1) + &
                            derivative_s(i)*weight_theta_phi*dbmod_s
                        second_value(2) = second_value(2) + weight_s(i)* &
                            derivative_theta(j)*weight_phi(k)*dbmod_s
                        second_value(3) = second_value(3) + weight_s(i)* &
                            weight_theta(j)*derivative_phi(k)*dbmod_s
                    end if
                    if (mode_secders == 2) then
                        second_value(4) = second_value(4) + weight_s(i)* &
                            second_theta(j)*weight_phi(k)*bmod
                        second_value(5) = second_value(5) + weight_s(i)* &
                            derivative_theta(j)*derivative_phi(k)*bmod
                        second_value(6) = second_value(6) + weight_s(i)* &
                            weight_theta(j)*second_phi(k)*bmod
                    end if
                end do
            end do
        end do
    end subroutine evaluate_symplectic_field_table

    subroutine splint_boozer_rk_table_device(s, theta, phi, aphi, daphi, &
            btheta, dbtheta, bphi, dbphi, bmod, dbmod)
        !$acc routine seq
        real(dp), intent(in) :: s, theta, phi
        real(dp), intent(out) :: aphi, daphi
        real(dp), intent(out) :: btheta, dbtheta, bphi, dbphi
        real(dp), intent(out) :: bmod, dbmod(3)

        real(sp) :: field_value(4), profile_value(6)

        call evaluate_rk_field_table(abs(s), theta, phi, field_value)
        call evaluate_rk_profile_table(abs(s), profile_value)
        bmod = real(field_value(1), dp)
        dbmod = real(field_value(2:4), dp)
        aphi = real(profile_value(1), dp)
        daphi = real(profile_value(2), dp)
        btheta = real(profile_value(3), dp)
        dbtheta = real(profile_value(4), dp)
        bphi = real(profile_value(5), dp)
        dbphi = real(profile_value(6), dp)
    end subroutine splint_boozer_rk_table_device

    subroutine splint_boozer_symplectic_table_device(s, theta, phi, &
            mode_secders, aphi, daphi, d2aphi, btheta, dbtheta, d2btheta, &
            bphi, dbphi, d2bphi, bmod, dbmod, d2bmod)
        !$acc routine seq
        real(dp), intent(in) :: s, theta, phi
        integer, intent(in) :: mode_secders
        real(dp), intent(out) :: aphi, daphi, d2aphi
        real(dp), intent(out) :: btheta, dbtheta, d2btheta
        real(dp), intent(out) :: bphi, dbphi, d2bphi, bmod
        real(dp), intent(out) :: dbmod(3), d2bmod(6)

        real(sp) :: field_value(4), field_second(6)
        real(sp) :: profile_value(6), profile_second(3)

        call evaluate_symplectic_field_table(abs(s), theta, phi, mode_secders, &
            field_value, field_second)
        call evaluate_symplectic_profile_table(abs(s), mode_secders, &
            profile_value, profile_second)
        bmod = real(field_value(1), dp)
        dbmod = real(field_value(2:4), dp)
        d2bmod = real(field_second, dp)
        aphi = real(profile_value(1), dp)
        daphi = real(profile_value(2), dp)
        d2aphi = real(profile_second(1), dp)
        btheta = real(profile_value(3), dp)
        dbtheta = real(profile_value(4), dp)
        d2btheta = real(profile_second(2), dp)
        bphi = real(profile_value(5), dp)
        dbphi = real(profile_value(6), dp)
        d2bphi = real(profile_second(3), dp)
    end subroutine splint_boozer_symplectic_table_device
end module boozer_rk_tables
