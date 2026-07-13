module jorek_model303_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_bezier, only: evaluate_jorek_geometry, locate_jorek_element
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
            a_r_phi_z, b_r_z_phi, element, st, ierr)
        type(jorek_restart_t), intent(in) :: data
        real(dp), intent(in) :: target_rz(2), phi
        real(dp), intent(out) :: a_r_phi_z(3), b_r_z_phi(3)
        integer, intent(out) :: element
        real(dp), intent(out) :: st(2)
        integer, intent(out) :: ierr

        a_r_phi_z = 0.0_dp
        b_r_z_phi = 0.0_dp
        element = 0
        st = 0.0_dp
        call locate_jorek_element(data, target_rz, element, st(1), st(2), ierr)
        if (ierr /= 0) return
        call evaluate_jorek_model303_a(data, element, st(1), st(2), phi, &
            a_r_phi_z, ierr)
        if (ierr /= 0) return
        call evaluate_jorek_model303_b(data, element, st(1), st(2), phi, &
            b_r_z_phi, ierr)
    end subroutine evaluate_jorek_model303_at

end module jorek_model303_field
