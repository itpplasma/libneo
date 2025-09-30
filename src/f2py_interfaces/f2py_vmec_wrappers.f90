module f2py_vmec_wrappers
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spline_vmec_sub, only : splint_vmec_data, vmec_field
    use vmec_field_tools, only : vmec_field_cylindrical
    implicit none
contains

    subroutine splint_vmec_data_wrapper(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                                        R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
        real(dp), intent(in) :: s, theta, varphi
        real(dp), intent(out) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota
        real(dp), intent(out) :: R, Z, alam
        real(dp), intent(out) :: dR_ds, dR_dt, dR_dp
        real(dp), intent(out) :: dZ_ds, dZ_dt, dZ_dp
        real(dp), intent(out) :: dl_ds, dl_dt, dl_dp

        call splint_vmec_data(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                              R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
    end subroutine splint_vmec_data_wrapper

    subroutine vmec_field_wrapper(s, theta, varphi, A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                                  sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, &
                                  Bcovar_s, Bcovar_vartheta, Bcovar_varphi)
        real(dp), intent(in) :: s, theta, varphi
        real(dp), intent(out) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota
        real(dp), intent(out) :: sqg, alam, dl_ds, dl_dt, dl_dp
        real(dp), intent(out) :: Bctrvr_vartheta, Bctrvr_varphi
        real(dp), intent(out) :: Bcovar_s, Bcovar_vartheta, Bcovar_varphi

        call vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                        sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, &
                        Bcovar_s, Bcovar_vartheta, Bcovar_varphi)
    end subroutine vmec_field_wrapper


    subroutine vmec_field_cylindrical_wrapper(s, theta, varphi, BR, Bphi, BZ, Bmag)
        real(dp), intent(in) :: s, theta, varphi
        real(dp), intent(out) :: BR, Bphi, BZ, Bmag

        call vmec_field_cylindrical(s, theta, varphi, BR, Bphi, BZ, Bmag)
    end subroutine vmec_field_cylindrical_wrapper

end module f2py_vmec_wrappers
