module vmec_field_tools

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spline_vmec_sub, only: splint_vmec_data, vmec_field

    implicit none

    real(dp), parameter :: cm_to_m = 1.0d-2
    real(dp), parameter :: gauss_to_tesla = 1.0d-4

contains

    subroutine vmec_field_cylindrical(s, theta, varphi, BR, Bphi, BZ, Bmag)
        real(dp), intent(in) :: s
        real(dp), intent(in) :: theta
        real(dp), intent(in) :: varphi
        real(dp), intent(out) :: BR
        real(dp), intent(out) :: Bphi
        real(dp), intent(out) :: BZ
        real(dp), intent(out) :: Bmag

        real(dp) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota
        real(dp) :: R_cm, Z_cm, alam
        real(dp) :: dR_ds_cm, dR_dt_cm, dR_dp_cm
        real(dp) :: dZ_ds_cm, dZ_dt_cm, dZ_dp_cm
        real(dp) :: dl_ds, dl_dt, dl_dp
        real(dp) :: sqg, Bctr_vartheta, Bctr_varphi
        real(dp) :: Bcov_s, Bcov_vartheta, Bcov_varphi
        real(dp) :: R, cosphi, sinphi
        real(dp) :: dR_ds, dR_dt, dR_dp
        real(dp) :: dZ_ds, dZ_dt, dZ_dp
        real(dp) :: cjac
        real(dp), dimension(3) :: e_s, e_theta, e_phi
        real(dp), dimension(3) :: e_vartheta, e_varphi
        real(dp), dimension(3) :: Bvec
        real(dp) :: Btotal

        call splint_vmec_data( s, theta, varphi, &
            A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
            R_cm, Z_cm, alam, &
            dR_ds_cm, dR_dt_cm, dR_dp_cm, &
            dZ_ds_cm, dZ_dt_cm, dZ_dp_cm, &
            dl_ds, dl_dt, dl_dp )

        call vmec_field( s, theta, varphi, &
            A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
            sqg, alam, dl_ds, dl_dt, dl_dp, &
            Bctr_vartheta, Bctr_varphi, &
            Bcov_s, Bcov_vartheta, Bcov_varphi )

        R = R_cm * cm_to_m
        cosphi = cos(varphi)
        sinphi = sin(varphi)

        dR_ds = dR_ds_cm * cm_to_m
        dR_dt = dR_dt_cm * cm_to_m
        dR_dp = dR_dp_cm * cm_to_m
        dZ_ds = dZ_ds_cm * cm_to_m
        dZ_dt = dZ_dt_cm * cm_to_m
        dZ_dp = dZ_dp_cm * cm_to_m

        e_s = [ dR_ds * cosphi, dR_ds * sinphi, dZ_ds ]
        e_theta = [ dR_dt * cosphi, dR_dt * sinphi, dZ_dt ]
        e_phi = [ dR_dp * cosphi - R * sinphi, &
                  dR_dp * sinphi + R * cosphi, &
                  dZ_dp ]

        cjac = 1.0_dp / (1.0_dp + dl_dt)
        e_vartheta = cjac * ( e_theta - dl_ds * e_s - dl_dp * e_phi )
        e_varphi = e_phi - dl_dp * cjac * e_theta

        Bvec = gauss_to_tesla * ( Bctr_vartheta * e_vartheta + Bctr_varphi * e_varphi )

        BR = Bvec(1) * cosphi + Bvec(2) * sinphi
        Bphi = -Bvec(1) * sinphi + Bvec(2) * cosphi
        BZ = Bvec(3)

        Btotal = sqrt(Bvec(1)**2 + Bvec(2)**2 + Bvec(3)**2)
        Bmag = Btotal
    end subroutine vmec_field_cylindrical

end module vmec_field_tools
