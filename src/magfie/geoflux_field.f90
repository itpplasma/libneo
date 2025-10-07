module geoflux_field

   use, intrinsic :: iso_fortran_env, only: dp => real64
   use geoflux_coordinates, only: init_geoflux_coordinates, geoflux_to_cyl, geoflux_get_axis
   use field_sub, only: field_eq, psif
   use field_eq_mod, only: psi_axis, reset_field_eq_state
   use input_files, only: gfile

   implicit none

   real(dp), save :: axis_R = 0.0_dp
   real(dp), save :: axis_Z = 0.0_dp
   logical, save :: initialized = .false.

contains

   subroutine spline_geoflux_data(geqdsk_file, ns_cache, ntheta_cache)
      character(len=*), intent(in) :: geqdsk_file
      integer, intent(in), optional :: ns_cache, ntheta_cache

      call reset_field_eq_state()
      gfile = trim(geqdsk_file)

      call init_geoflux_coordinates(geqdsk_file, ns_cache, ntheta_cache)
      call geoflux_get_axis(axis_R, axis_Z)

      initialized = .true.
   end subroutine spline_geoflux_data

   subroutine splint_geoflux_field(s_tor, theta_geo, phi, Acov, hcov, Bmod, sqgBctr)
      real(dp), intent(in) :: s_tor, theta_geo, phi
      real(dp), intent(out) :: Acov(3), hcov(3), Bmod
      real(dp), intent(out), optional :: sqgBctr(3)

      real(dp) :: x_geo(3)
      real(dp) :: x_cyl(3)
      real(dp) :: jac(3, 3)
      real(dp) :: dRds, dRdt, dZds, dZdt
      real(dp) :: det_s_theta, sqrtg
      real(dp) :: Br, Bphi, Bz
      real(dp) :: dBrdR, dBrdp, dBrdZ
      real(dp) :: dBpdR, dBpdp, dBpdZ
      real(dp) :: dBzdR, dBzdp, dBzdZ
      real(dp) :: Bcov_s, Bcov_theta, Bcov_phi
      real(dp) :: psi_local

      if (.not. initialized) then
         error stop "geoflux_field: call spline_geoflux_data before splint_geoflux_field"
      end if

      x_geo = [s_tor, theta_geo, phi]
      call geoflux_to_cyl(x_geo, x_cyl, jac)

      call field_eq(x_cyl(1), x_cyl(2), x_cyl(3), Br, Bphi, Bz, &
                    dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

      psi_local = psif

      Bmod = sqrt(Br*Br + Bphi*Bphi + Bz*Bz)

      dRds = jac(1, 1)
      dRdt = jac(1, 2)
      dZds = jac(3, 1)
      dZdt = jac(3, 2)

      det_s_theta = dRds*dZdt - dRdt*dZds
      sqrtg = abs(x_cyl(1)*det_s_theta)

      Bcov_s = Br*dRds + Bz*dZds
      Bcov_theta = Br*dRdt + Bz*dZdt
      Bcov_phi = Bphi*x_cyl(1)

      if (Bmod > 0.0_dp) then
         hcov(1) = Bcov_s/Bmod
         hcov(2) = Bcov_theta/Bmod
         hcov(3) = Bcov_phi/Bmod
      else
         hcov = 0.0_dp
      end if

      Acov(1) = 0.0_dp
      Acov(2) = 0.0_dp
      Acov(3) = psi_local - psi_axis

      if (present(sqgBctr)) then
         sqgBctr = 0.0_dp
         sqgBctr(1) = sqrtg
      end if
   end subroutine splint_geoflux_field

end module geoflux_field
