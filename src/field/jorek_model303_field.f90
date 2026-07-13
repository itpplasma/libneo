module jorek_model303_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: jorek_locator_t, evaluate_jorek_geometry, &
        locate_jorek_element, locate_jorek_element_indexed, &
        is_jorek_axis_target
    use jorek_field_values, only: evaluate_jorek_variable
    use jorek_restart, only: jorek_restart_t

    implicit none
    private

    public :: evaluate_jorek_model303_a, evaluate_jorek_model303_b, &
        evaluate_jorek_model303_at

contains

    pure subroutine evaluate_jorek_model303_a(data, element, s, t, phi, &
            a_r_phi_z, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element
        real(dp), intent(in) :: s, t, phi
        real(dp), intent(out) :: a_r_phi_z(3)
        integer, intent(out) :: ierr

        real(dp) :: rz(2), rz_st(2, 2), psi, psi_st_phi(3)

        a_r_phi_z = 0.0_dp
        ierr = 0
        if (data%jorek_model /= 303) then
            ierr = 7
            return
        end if
        call evaluate_jorek_geometry(data, element, s, t, rz, rz_st, ierr)
        if (ierr /= 0) return
        call evaluate_jorek_variable(data, element, 1, s, t, phi, psi, &
            psi_st_phi, ierr)
        if (ierr /= 0) return
        if (rz(1) <= 0.0_dp) then
            ierr = 8
            return
        end if

        a_r_phi_z = [0.0_dp, -psi, -data%F0*log(rz(1))]
    end subroutine evaluate_jorek_model303_a

    pure subroutine evaluate_jorek_model303_b(data, element, s, t, phi, &
            b_r_z_phi, ierr)
        type(jorek_restart_t), intent(in) :: data
        integer, intent(in) :: element
        real(dp), intent(in) :: s, t, phi
        real(dp), intent(out) :: b_r_z_phi(3)
        integer, intent(out) :: ierr

        real(dp) :: rz(2), rz_st(2, 2), psi, psi_st_phi(3)
        real(dp) :: jacobian, jacobian_scale, psi_r, psi_z

        b_r_z_phi = 0.0_dp
        ierr = 0
        if (data%jorek_model /= 303) then
            ierr = 7
            return
        end if
        call evaluate_jorek_geometry(data, element, s, t, rz, rz_st, ierr)
        if (ierr /= 0) return
        call evaluate_jorek_variable(data, element, 1, s, t, phi, psi, &
            psi_st_phi, ierr)
        if (ierr /= 0) return

        jacobian = rz_st(1, 1)*rz_st(2, 2) - rz_st(1, 2)*rz_st(2, 1)
        jacobian_scale = max(abs(rz_st(1, 1)*rz_st(2, 2)), &
            abs(rz_st(1, 2)*rz_st(2, 1)), tiny(1.0_dp))
        if (rz(1) <= 0.0_dp &
            .or. abs(jacobian) <= 16.0_dp*epsilon(1.0_dp)*jacobian_scale) then
            ierr = 8
            return
        end if

        psi_r = (psi_st_phi(1)*rz_st(2, 2) &
            - psi_st_phi(2)*rz_st(2, 1))/jacobian
        psi_z = (-psi_st_phi(1)*rz_st(1, 2) &
            + psi_st_phi(2)*rz_st(1, 1))/jacobian
        b_r_z_phi = [psi_z, -psi_r, data%F0]/rz(1)
    end subroutine evaluate_jorek_model303_b

    pure subroutine evaluate_jorek_model303_at(data, target_rz, phi, &
            a_r_phi_z, b_r_z_phi, element, st, ierr, locator)
        type(jorek_restart_t), intent(in) :: data
        real(dp), intent(in) :: target_rz(2), phi
        real(dp), intent(out) :: a_r_phi_z(3), b_r_z_phi(3)
        integer, intent(out) :: element
        real(dp), intent(out) :: st(2)
        integer, intent(out) :: ierr
        type(jorek_locator_t), intent(in), optional :: locator

        a_r_phi_z = 0.0_dp
        b_r_z_phi = 0.0_dp
        element = 0
        st = 0.0_dp
        if (present(locator)) then
            if (is_jorek_axis_target(locator, target_rz)) then
                call evaluate_jorek_model303_axis(data, locator, phi, &
                    a_r_phi_z, b_r_z_phi, element, st, ierr)
                return
            end if
            call locate_jorek_element_indexed(data, locator, target_rz, &
                element, st(1), st(2), ierr)
        else
            call locate_jorek_element(data, target_rz, element, st(1), st(2), &
                ierr)
        end if
        if (ierr /= 0) return
        call evaluate_jorek_model303_a(data, element, st(1), st(2), phi, &
            a_r_phi_z, ierr)
        if (ierr /= 0) return
        call evaluate_jorek_model303_b(data, element, st(1), st(2), phi, &
            b_r_z_phi, ierr)
    end subroutine evaluate_jorek_model303_at

    pure subroutine evaluate_jorek_model303_axis(data, locator, phi, &
            a_r_phi_z, b_r_z_phi, element, st, ierr)
        type(jorek_restart_t), intent(in) :: data
        type(jorek_locator_t), intent(in) :: locator
        real(dp), intent(in) :: phi
        real(dp), intent(out) :: a_r_phi_z(3), b_r_z_phi(3)
        integer, intent(out) :: element
        real(dp), intent(out) :: st(2)
        integer, intent(out) :: ierr

        real(dp), parameter :: delta(2) = [1.0e-2_dp, 1.0e-3_dp]
        real(dp), allocatable :: a_sample(:, :, :), angles(:), b_sample(:, :, :)
        real(dp), allocatable :: weights(:)
        real(dp) :: a_mean(3, 2), b_mean(3, 2), rz(2), rz_st(2, 2)
        real(dp) :: weight_sum
        integer :: i, level, n_axis

        a_r_phi_z = 0.0_dp
        b_r_z_phi = 0.0_dp
        element = 0
        st = 0.0_dp
        ierr = 8
        if (data%jorek_model /= 303) then
            ierr = 7
            return
        end if
        if (.not. locator%has_axis .or. locator%axis_rz(1) <= 0.0_dp) return
        if (.not. allocated(locator%axis_elements) &
            .or. .not. allocated(locator%axis_sides)) return
        n_axis = size(locator%axis_elements)
        if (n_axis < 1 .or. size(locator%axis_sides) /= n_axis) return
        allocate (a_sample(3, n_axis, 2), b_sample(3, n_axis, 2))
        allocate (angles(n_axis), weights(n_axis))
        do level = 1, 2
            do i = 1, n_axis
                call inward_axis_coordinates(locator%axis_sides(i), &
                    delta(level), st)
                call evaluate_jorek_geometry(data, locator%axis_elements(i), &
                    st(1), st(2), rz, rz_st, ierr)
                if (ierr /= 0) return
                if (level == 2) angles(i) = atan2(rz(2) - locator%axis_rz(2), &
                    rz(1) - locator%axis_rz(1))
                call evaluate_jorek_model303_a(data, &
                    locator%axis_elements(i), st(1), st(2), phi, &
                    a_sample(:, i, level), ierr)
                if (ierr /= 0) return
                call evaluate_jorek_model303_b(data, &
                    locator%axis_elements(i), st(1), st(2), phi, &
                    b_sample(:, i, level), ierr)
                if (ierr /= 0) return
                b_sample(1:2, i, level) = &
                    rz(1)*b_sample(1:2, i, level)
            end do
        end do
        call angular_voronoi_weights(angles, weights)
        weight_sum = sum(weights)
        if (weight_sum <= 0.0_dp) then
            ierr = 8
            return
        end if
        do level = 1, 2
            do i = 1, 3
                a_mean(i, level) = sum(weights*a_sample(i, :, level))/weight_sum
                b_mean(i, level) = sum(weights*b_sample(i, :, level))/weight_sum
            end do
        end do
        a_r_phi_z = (10.0_dp*a_mean(:, 2) - a_mean(:, 1))/9.0_dp
        b_r_z_phi = (10.0_dp*b_mean(:, 2) - b_mean(:, 1))/9.0_dp
        a_r_phi_z(1) = 0.0_dp
        a_r_phi_z(3) = -data%F0*log(locator%axis_rz(1))
        b_r_z_phi(1:2) = b_r_z_phi(1:2)/locator%axis_rz(1)
        b_r_z_phi(3) = data%F0/locator%axis_rz(1)
        element = locator%axis_elements(1)
        call inward_axis_coordinates(locator%axis_sides(1), 0.0_dp, st)
        ierr = 0
    end subroutine evaluate_jorek_model303_axis

    pure subroutine inward_axis_coordinates(side, delta, st)
        integer, intent(in) :: side
        real(dp), intent(in) :: delta
        real(dp), intent(out) :: st(2)

        select case (side)
        case (1)
            st = [delta, 0.5_dp]
        case (2)
            st = [1.0_dp - delta, 0.5_dp]
        case (3)
            st = [0.5_dp, delta]
        case (4)
            st = [0.5_dp, 1.0_dp - delta]
        case default
            st = 0.0_dp
        end select
    end subroutine inward_axis_coordinates

    pure subroutine angular_voronoi_weights(angles, weights)
        real(dp), intent(in) :: angles(:)
        real(dp), intent(out) :: weights(size(angles))

        real(dp), parameter :: two_pi = 2.0_dp*acos(-1.0_dp)
        real(dp) :: backward, difference, forward
        integer :: i, j

        if (size(angles) == 1) then
            weights = 1.0_dp
            return
        end if
        do i = 1, size(angles)
            forward = two_pi
            backward = two_pi
            do j = 1, size(angles)
                if (j == i) cycle
                difference = modulo(angles(j) - angles(i), two_pi)
                if (difference > 0.0_dp) forward = min(forward, difference)
                difference = modulo(angles(i) - angles(j), two_pi)
                if (difference > 0.0_dp) backward = min(backward, difference)
            end do
            weights(i) = 0.5_dp*(forward + backward)
        end do
    end subroutine angular_voronoi_weights

end module jorek_model303_field
