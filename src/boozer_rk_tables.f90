module boozer_rk_tables
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    public :: rk_field_table, rk_profile_table
    public :: rk_num_points, rk_x_min, rk_h_step, rk_inv_h_step, rk_period
    public :: rk_tables_ready
    public :: splint_boozer_rk_table_device

    ! Keep the compact Landreman-compatible tables in a small module of their
    ! own. NVHPC 26.5 otherwise combines these declarations with the generic
    ! spline descriptors and emits an invalid device-side static initializer.
    real(dp), allocatable, save :: rk_field_table(:)
    real(dp), allocatable, save :: rk_profile_table(:)
    integer, save :: rk_num_points(3)
    real(dp), save :: rk_x_min(3)
    real(dp), save :: rk_h_step(3)
    real(dp), save :: rk_inv_h_step(3)
    real(dp), save :: rk_period(3)
    logical, save :: rk_tables_ready = .false.
    !$acc declare create(rk_field_table, rk_profile_table, rk_num_points, &
    !$acc& rk_x_min, rk_h_step, rk_inv_h_step, rk_period, rk_tables_ready)

contains

    subroutine cubic_table_location(x, idim, first, weight)
        !$acc routine seq
        real(dp), intent(in) :: x
        integer, intent(in) :: idim
        integer, intent(out) :: first
        real(dp), intent(out) :: weight(4)

        real(dp) :: x_eval, x_grid, relative
        integer :: first_zero

        x_eval = x
        if (idim > 1) then
            x_eval = rk_x_min(idim) + modulo( &
                x - rk_x_min(idim), rk_period(idim))
        end if
        x_grid = (x_eval - rk_x_min(idim))*rk_inv_h_step(idim)
        first_zero = 3*(int(x_grid)/3)
        first_zero = max(0, min(first_zero, rk_num_points(idim) - 4))
        relative = x_grid - real(first_zero, dp)
        first = first_zero + 1

        weight(1) = (1.0_dp - relative)*(2.0_dp - relative)* &
            (3.0_dp - relative)/6.0_dp
        weight(2) = relative*(2.0_dp - relative)* &
            (3.0_dp - relative)/2.0_dp
        weight(3) = relative*(relative - 1.0_dp)* &
            (3.0_dp - relative)/2.0_dp
        weight(4) = relative*(relative - 1.0_dp)* &
            (relative - 2.0_dp)/6.0_dp
    end subroutine cubic_table_location

    subroutine evaluate_rk_profile_table(s, value)
        !$acc routine seq
        real(dp), intent(in) :: s
        real(dp), intent(out) :: value(6)

        integer :: first, i, iq, table_index
        real(dp) :: weight(4)

        call cubic_table_location(s, 1, first, weight)
        value = 0.0_dp
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
        real(dp), intent(out) :: value(4)

        integer :: first_s, first_theta, first_phi
        integer :: i, j, k, iq, table_index
        real(dp) :: weight_s(4), weight_theta(4), weight_phi(4), weight

        call cubic_table_location(s, 1, first_s, weight_s)
        call cubic_table_location(theta, 2, first_theta, weight_theta)
        call cubic_table_location(phi, 3, first_phi, weight_phi)
        value = 0.0_dp
        do k = 1, 4
            do j = 1, 4
                do i = 1, 4
                    weight = weight_s(i)*weight_theta(j)*weight_phi(k)
                    do iq = 1, 4
                        table_index = iq + 4*((first_s + i - 2) + &
                            rk_num_points(1)*((first_theta + j - 2) + &
                            rk_num_points(2)*(first_phi + k - 2)))
                        value(iq) = value(iq) + &
                            weight*rk_field_table(table_index)
                    end do
                end do
            end do
        end do
    end subroutine evaluate_rk_field_table

    subroutine splint_boozer_rk_table_device(s, theta, phi, aphi, daphi, &
            btheta, dbtheta, bphi, dbphi, bmod, dbmod)
        !$acc routine seq
        real(dp), intent(in) :: s, theta, phi
        real(dp), intent(out) :: aphi, daphi
        real(dp), intent(out) :: btheta, dbtheta, bphi, dbphi
        real(dp), intent(out) :: bmod, dbmod(3)

        real(dp) :: field_value(4), profile_value(6)

        call evaluate_rk_field_table(abs(s), theta, phi, field_value)
        call evaluate_rk_profile_table(abs(s), profile_value)
        bmod = field_value(1)
        dbmod = field_value(2:4)
        aphi = profile_value(1)
        daphi = profile_value(2)
        btheta = profile_value(3)
        dbtheta = profile_value(4)
        bphi = profile_value(5)
        dbphi = profile_value(6)
    end subroutine splint_boozer_rk_table_device
end module boozer_rk_tables
