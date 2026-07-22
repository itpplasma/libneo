module spline_vmec_sub
   use iso_fortran_env, only: dp => real64
   use spl_three_to_five_sub

   implicit none

contains

   subroutine spline_vmec_data
      use vmec_alloc_sub, only: new_allocate_vmec_stuff, new_deallocate_vmec_stuff
      call new_allocate_vmec_stuff
      call initialize_vmec_data_and_arrays
      call setup_poloidal_flux_splines
      call setup_angular_grid_and_fourier_synthesis
      call new_deallocate_vmec_stuff

      call spline_over_phi
      call spline_over_theta
      call spline_over_s
   end subroutine spline_vmec_data

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine initialize_vmec_data_and_arrays
      use new_vmec_stuff_mod, only: ns_A, ns_s, ns_tp, rmnc, zmns, almns, rmns, zmnc, almnc, &
                                    aiota, phi, sps, axm, axn, s, nsurfm, nstrm, kpar
      use vmecin_sub, only: vmecin
      use vector_potentail_mod, only: ns, hs, torflux

      print *, 'Splining VMEC data: ns_A = ', ns_A, '  ns_s = ', ns_s, '  ns_tp = ', ns_tp

      call vmecin(rmnc, zmns, almns, rmns, zmnc, almnc, aiota, phi, sps, axm, axn, s, &
                  nsurfm, nstrm, kpar, torflux)

      ns = kpar + 1
      hs = s(2) - s(1)
   end subroutine initialize_vmec_data_and_arrays

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine setup_poloidal_flux_splines
      use new_vmec_stuff_mod, only: ns_A, aiota
      use vector_potentail_mod, only: ns, hs, torflux, sA_phi
      use spl_three_to_five_sub, only: spl_reg

      integer :: i, is, k
      real(dp), dimension(:, :), allocatable :: splcoe

      allocate (splcoe(0:ns_A, ns))

      splcoe(0, :) = aiota

      call spl_reg(ns_A - 1, ns, hs, splcoe(0:ns_A - 1, :))

      do i = ns_A, 1, -1
         splcoe(i, :) = splcoe(i - 1, :)/dble(i)
      end do

      splcoe(0, 1) = 0.d0
      do is = 1, ns - 1
         splcoe(0, is + 1) = splcoe(ns_A, is)
         do k = ns_A - 1, 0, -1
            splcoe(0, is + 1) = splcoe(k, is) + hs*splcoe(0, is + 1)
         end do
      end do

      ! A_phi = -chi = -(poloidal flux/2pi) = -torflux * integral(iota ds), the
      ! toroidal covariant component A_zeta. The toroidal flux psi_tor sits on the
      ! poloidal component instead: A_theta = torflux*s (interpolate_vector_potential).
      ! For VMEC input, sign(torflux)=sign(signgs*phi_edge).
      if (allocated(sA_phi)) deallocate (sA_phi)
      allocate (sA_phi(ns_A + 1, ns))
      do k = 0, ns_A
         sA_phi(k + 1, :) = -torflux*splcoe(k, :)
      end do

      deallocate (splcoe)
   end subroutine setup_poloidal_flux_splines

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine setup_angular_grid_and_fourier_synthesis
      use new_vmec_stuff_mod, only: axm, axn, multharm, nper, nstrm, &
                                    n_theta, n_phi, h_theta, h_phi, sR, sZ, slam, ns_s, ns_tp
      use vector_potentail_mod, only: ns
      use vmec_alloc_sub, only: new_deallocate_vmec_stuff

      real(dp), dimension(:, :), allocatable :: almnc_rho, rmnc_rho, zmnc_rho
      real(dp), dimension(:, :), allocatable :: almns_rho, rmns_rho, zmns_rho

      integer :: i, m, n, is, i_theta, i_phi, m_max, n_max
      integer :: nsize_exp_imt, nsize_exp_inp, iexpt, iexpp
      real(dp) :: twopi, cosphase, sinphase
      complex(8) :: base_exp_imt, base_exp_inp, base_exp_inp_inv, expphase
      complex(8), dimension(:), allocatable :: exp_imt, exp_inp

      call perform_axis_healing(almnc_rho, rmnc_rho, zmnc_rho, almns_rho, rmns_rho, zmns_rho)

      ! Setup angular grid
      m_max = nint(maxval(axm))
      n_max = nint(maxval(axn))

      print *, 'VMEC ns = ', ns, ' m_max = ', m_max, ' n_max = ', n_max

      n_theta = m_max*multharm + 1
      n_phi = n_max*multharm + 1
      twopi = 8.d0*atan2(1.d0, 1.d0)
      h_theta = twopi/dble(n_theta - 1)
      h_phi = twopi/dble((n_phi - 1)*nper)

      ! Setup exponential arrays
      nsize_exp_imt = (n_theta - 1)*m_max
      nsize_exp_inp = (n_phi - 1)*n_max

      allocate (exp_imt(0:nsize_exp_imt), exp_inp(-nsize_exp_inp:nsize_exp_inp))

      base_exp_imt = exp(cmplx(0.d0, h_theta, kind=kind(0d0)))
      base_exp_inp = exp(cmplx(0.d0, h_phi, kind=kind(0d0)))
      base_exp_inp_inv = (1.d0, 0.d0)/base_exp_inp
      exp_imt(0) = (1.d0, 0.d0)
      exp_inp(0) = (1.d0, 0.d0)

      do i = 1, nsize_exp_imt
         exp_imt(i) = exp_imt(i - 1)*base_exp_imt
      end do

      do i = 1, nsize_exp_inp
         exp_inp(i) = exp_inp(i - 1)*base_exp_inp
         exp_inp(-i) = exp_inp(1 - i)*base_exp_inp_inv
      end do

      ! Allocate spline arrays
      if (allocated(sR)) deallocate (sR)
      allocate (sR(ns_s + 1, ns_tp + 1, ns_tp + 1, ns, n_theta, n_phi))
      if (allocated(sZ)) deallocate (sZ)
      allocate (sZ(ns_s + 1, ns_tp + 1, ns_tp + 1, ns, n_theta, n_phi))
      if (allocated(slam)) deallocate (slam)
      allocate (slam(ns_s + 1, ns_tp + 1, ns_tp + 1, ns, n_theta, n_phi))

      sR(1, 1, 1, :, :, :) = 0.d0
      sZ(1, 1, 1, :, :, :) = 0.d0
      slam(1, 1, 1, :, :, :) = 0.d0

      ! Fourier synthesis to real space
!$omp parallel private(m, n, i_theta, i_phi, i, is, iexpt, iexpp, expphase, cosphase, sinphase)
!$omp do
      do i_theta = 1, n_theta
         do i_phi = 1, n_phi
            do i = 1, nstrm
               m = nint(axm(i))
               n = nint(axn(i))
               iexpt = m*(i_theta - 1)
               iexpp = n*(i_phi - 1)
               expphase = exp_imt(iexpt)*exp_inp(-iexpp)
               cosphase = dble(expphase)
               sinphase = aimag(expphase)
               do is = 1, ns
                  sR(1, 1, 1, is, i_theta, i_phi) = sR(1, 1, 1, is, i_theta, i_phi) &
                                                    + rmnc_rho(i, is - 1)*cosphase &
                                                    + rmns_rho(i, is - 1)*sinphase
                  sZ(1, 1, 1, is, i_theta, i_phi) = sZ(1, 1, 1, is, i_theta, i_phi) &
                                                    + zmnc_rho(i, is - 1)*cosphase &
                                                    + zmns_rho(i, is - 1)*sinphase
                  slam(1, 1, 1, is, i_theta, i_phi) = slam(1, 1, 1, is, i_theta, i_phi) &
                                                      + almnc_rho(i, is - 1)*cosphase &
                                                      + almns_rho(i, is - 1)*sinphase
               end do
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      deallocate (exp_imt, exp_inp)
      deallocate (almnc_rho, rmnc_rho, zmnc_rho, almns_rho, rmns_rho, zmns_rho)
   end subroutine setup_angular_grid_and_fourier_synthesis

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine perform_axis_healing(almnc_rho, rmnc_rho, zmnc_rho, almns_rho, rmns_rho, zmns_rho)
      use new_vmec_stuff_mod, only: rmnc, zmns, almns, rmns, zmnc, almnc, &
                                    axm, axn, nstrm, nper, axis_healing_polyfit_degree, &
                                    raxis_cc, zaxis_cs, raxis_cs, zaxis_cc, ntor_axis
      use vector_potentail_mod, only: ns

      real(dp), dimension(:, :), allocatable, intent(out) :: almnc_rho, rmnc_rho, zmnc_rho
      real(dp), dimension(:, :), allocatable, intent(out) :: almns_rho, rmns_rho, zmns_rho
      integer :: i, m, nrho, nheal, i_anchor, n_idx
      real(dp) :: anc_rc, anc_zs, anc_rs, anc_zc
      character(len=16) :: mode, boundary

      nrho = ns
      allocate (almnc_rho(nstrm, 0:nrho - 1), rmnc_rho(nstrm, 0:nrho - 1), zmnc_rho(nstrm, 0:nrho - 1))
      allocate (almns_rho(nstrm, 0:nrho - 1), rmns_rho(nstrm, 0:nrho - 1), zmns_rho(nstrm, 0:nrho - 1))

      call resolve_axis_healing_mode(mode, boundary)

      select case (trim(mode))
      case ('polyfit')
         i_anchor = axis_anchor_index(ns)
         print *, 'VMEC axis healing: polyfit, rho**m * poly(s) below s = ', &
            dble(i_anchor - 1)/dble(ns - 1), ' degree ', axis_healing_polyfit_degree
         do i = 1, nstrm
            m = nint(abs(axm(i)))
            ! Pin the m=0 geometry amplitudes to the exact magnetic axis
            ! (raxis_cc/zaxis_cs); lambda has no axis value, m/=0 ignores it.
            anc_rc = 0.0d0; anc_zs = 0.0d0; anc_rs = 0.0d0; anc_zc = 0.0d0
            if (m == 0 .and. allocated(raxis_cc)) then
               n_idx = nint(axn(i))/nper
               if (n_idx >= 0 .and. n_idx <= ntor_axis) then
                  anc_rc = raxis_cc(n_idx); anc_zs = zaxis_cs(n_idx)
                  anc_rs = raxis_cs(n_idx); anc_zc = zaxis_cc(n_idx)
               end if
            end if
            call s_to_rho_polyfit(m, ns, nrho, i_anchor, axis_healing_polyfit_degree, rmnc(i, :), rmnc_rho(i, :), anchor=anc_rc)
            call s_to_rho_polyfit(m, ns, nrho, i_anchor, axis_healing_polyfit_degree, zmnc(i, :), zmnc_rho(i, :), anchor=anc_zc)
            call s_to_rho_polyfit(m, ns, nrho, i_anchor, axis_healing_polyfit_degree, almnc(i, :), almnc_rho(i, :))
            call s_to_rho_polyfit(m, ns, nrho, i_anchor, axis_healing_polyfit_degree, rmns(i, :), rmns_rho(i, :), anchor=anc_rs)
            call s_to_rho_polyfit(m, ns, nrho, i_anchor, axis_healing_polyfit_degree, zmns(i, :), zmns_rho(i, :), anchor=anc_zs)
            call s_to_rho_polyfit(m, ns, nrho, i_anchor, axis_healing_polyfit_degree, almns(i, :), almns_rho(i, :))
         end do

      case ('powerlaw')
         i_anchor = axis_anchor_index(ns)
         print *, 'VMEC axis healing: powerlaw, rho**m below s = ', &
            dble(i_anchor - 1)/dble(ns - 1)
         do i = 1, nstrm
            m = nint(abs(axm(i)))
            call s_to_rho_power_law(m, ns, nrho, i_anchor, rmnc(i, :), rmnc_rho(i, :))
            call s_to_rho_power_law(m, ns, nrho, i_anchor, zmnc(i, :), zmnc_rho(i, :))
            call s_to_rho_power_law(m, ns, nrho, i_anchor, almnc(i, :), almnc_rho(i, :))
            call s_to_rho_power_law(m, ns, nrho, i_anchor, rmns(i, :), rmns_rho(i, :))
            call s_to_rho_power_law(m, ns, nrho, i_anchor, zmns(i, :), zmns_rho(i, :))
            call s_to_rho_power_law(m, ns, nrho, i_anchor, almns(i, :), almns_rho(i, :))
         end do

      case ('legacy', 'legacy_adaptive')
         do i = 1, nstrm
            m = nint(abs(axm(i)))
            if (trim(mode) == 'legacy') then
               nheal = min(m, 4)
            else
               call determine_nheal_for_axis(m, ns, rmnc(i, :), nheal)
            end if
            call s_to_rho_healaxis(m, ns, nrho, nheal, rmnc(i, :), rmnc_rho(i, :))
            call s_to_rho_healaxis(m, ns, nrho, nheal, zmnc(i, :), zmnc_rho(i, :))
            call s_to_rho_healaxis(m, ns, nrho, nheal, almnc(i, :), almnc_rho(i, :))
            call s_to_rho_healaxis(m, ns, nrho, nheal, rmns(i, :), rmns_rho(i, :))
            call s_to_rho_healaxis(m, ns, nrho, nheal, zmns(i, :), zmns_rho(i, :))
            call s_to_rho_healaxis(m, ns, nrho, nheal, almns(i, :), almns_rho(i, :))
         end do
      end select
   end subroutine perform_axis_healing

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine resolve_axis_healing_mode(mode, boundary)
      !> Resolves the axis-healing selector. The preferred inputs are the strings
      !> axis_healing ('legacy'|'legacy_adaptive'|'powerlaw'|'polyfit').
      !> When axis_healing is empty the deprecated logical switches are mapped to
      !> a mode and a deprecation warning is emitted. Unknown values stop the run.
      use new_vmec_stuff_mod, only: axis_healing, old_axis_healing_boundary, &
                                    axis_healing_power_law, axis_healing_polyfit

      character(len=*), intent(out) :: mode, boundary
      logical :: from_legacy

      from_legacy = (len_trim(axis_healing) == 0)
      if (from_legacy) then
         if (axis_healing_polyfit) then
            mode = 'polyfit'
         else if (axis_healing_power_law) then
            mode = 'powerlaw'
         else if (old_axis_healing_boundary) then
            mode = 'legacy'
         else
            mode = 'legacy_adaptive'
         end if
      else
         mode = to_lower(adjustl(axis_healing))
      end if

      select case (trim(mode))
      case ('legacy')
         boundary = 'fixed'
      case ('legacy_adaptive')
         boundary = 'adaptive'
      case ('powerlaw', 'polyfit')
         boundary = ''
      case default
         print *, 'ERROR: unknown axis_healing = ''', trim(axis_healing), ''''
         error stop 'unknown axis_healing mode (use legacy|legacy_adaptive|powerlaw|polyfit)'
      end select

      if (from_legacy) then
         print *, 'WARNING: axis_healing not set; mapped legacy switches to axis_healing=''', &
            trim(mode), '''. The old_axis_healing*/axis_healing_power_law/axis_healing_polyfit'
         print *, '         switches are deprecated; set axis_healing=''', trim(mode), ''' instead.'
         print *, '         These switches will become errors in the next SIMPLE release.'
      end if
      if (trim(mode) /= 'polyfit') then
         print *, 'NOTE: axis_healing=''', trim(mode), '''; ''polyfit'' is recommended and will'
         print *, '      become the default. Please test it.'
      end if
   end subroutine resolve_axis_healing_mode

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   pure function to_lower(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: t
      integer :: i, c
      t = s
      do i = 1, len_trim(s)
         c = iachar(s(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) t(i:i) = achar(c + 32)
      end do
   end function to_lower

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine spline_over_phi
      use new_vmec_stuff_mod, only: ns_tp, n_theta, n_phi, h_phi, sR, sZ, slam
      use vector_potentail_mod, only: ns
      integer :: is, i_theta

      if (n_phi == 1) then
         print *, 'Spline not supported for a Phi period of 1, exiting...'
         call exit(-1)
      end if

!$omp parallel private(is, i_theta)
!$omp do
      do is = 1, ns
         do i_theta = 1, n_theta
            call spline_1d_periodic(sR(1, 1, :, is, i_theta, :), ns_tp, h_phi)
            call spline_1d_periodic(sZ(1, 1, :, is, i_theta, :), ns_tp, h_phi)
            call spline_1d_periodic(slam(1, 1, :, is, i_theta, :), ns_tp, h_phi)
         end do
      end do
!$omp end do
!$omp end parallel
   end subroutine spline_over_phi

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine spline_over_theta
      use new_vmec_stuff_mod, only: ns_tp, n_phi, h_theta, sR, sZ, slam
      use vector_potentail_mod, only: ns
      integer :: is, i_phi, isp

!$omp parallel private(is, i_phi, isp)
!$omp do
      do is = 1, ns
         do i_phi = 1, n_phi
            do isp = 1, ns_tp + 1
               call spline_1d_periodic(sR(1, :, isp, is, :, i_phi), ns_tp, h_theta)
               call spline_1d_periodic(sZ(1, :, isp, is, :, i_phi), ns_tp, h_theta)
               call spline_1d_periodic(slam(1, :, isp, is, :, i_phi), ns_tp, h_theta)
            end do
         end do
      end do
!$omp end do
!$omp end parallel
   end subroutine spline_over_theta

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine spline_over_s
      use new_vmec_stuff_mod, only: ns_s, ns_tp, n_theta, n_phi, sR, sZ, slam
      use vector_potentail_mod, only: hs
      integer :: i_theta, i_phi, ist, isp

!$omp parallel private(i_theta, i_phi, ist, isp)
!$omp do
      do i_theta = 1, n_theta
         do i_phi = 1, n_phi
            do ist = 1, ns_tp + 1
               do isp = 1, ns_tp + 1
                  call spline_1d_regular(sR(:, ist, isp, :, i_theta, i_phi), ns_s, hs)
                  call spline_1d_regular(sZ(:, ist, isp, :, i_theta, i_phi), ns_s, hs)
                  call spline_1d_regular(slam(:, ist, isp, :, i_theta, i_phi), ns_s, hs)
               end do
            end do
         end do
      end do
!$omp end do
!$omp end parallel
   end subroutine spline_over_s

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine spline_1d_periodic(data_array, ns_order, h_step)
      use spl_three_to_five_sub, only: spl_per

      real(dp), dimension(:, :), intent(inout) :: data_array
      integer, intent(in) :: ns_order
      real(dp), intent(in) :: h_step

      integer :: n_points
      real(dp), dimension(0:ns_order, size(data_array, 2)) :: splcoe

      n_points = size(data_array, 2)
      splcoe(0, :) = data_array(1, :)
      call spl_per(ns_order, n_points, h_step, splcoe)
      data_array(1:ns_order + 1, :) = splcoe(0:ns_order, :)
   end subroutine spline_1d_periodic

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine spline_1d_regular(data_array, ns_order, h_step)
      use spl_three_to_five_sub, only: spl_reg

      real(dp), dimension(:, :), intent(inout) :: data_array
      integer, intent(in) :: ns_order
      real(dp), intent(in) :: h_step

      integer :: n_points
      real(dp), dimension(0:ns_order, size(data_array, 2)) :: splcoe

      n_points = size(data_array, 2)
      splcoe(0, :) = data_array(1, :)
      call spl_reg(ns_order, n_points, h_step, splcoe)
      data_array(1:ns_order + 1, :) = splcoe(0:ns_order, :)
   end subroutine spline_1d_regular

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine deallocate_vmec_spline(mode)

      use new_vmec_stuff_mod, only: sR, sZ, slam

      integer :: mode

      if (mode .eq. 0) then
         deallocate (sR, sZ, slam)
      elseif (mode .eq. 1) then
         deallocate (sR, sZ)
      elseif (mode .eq. 2) then
         deallocate (slam)
      else
         print *, 'deallocate_vmec_spline: unknown mode'
      end if

   end subroutine deallocate_vmec_spline

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine normalize_coordinates(s, theta, varphi, ds, is, dtheta, i_theta, dphi, i_phi)
      use new_vmec_stuff_mod, only: n_theta, n_phi, h_theta, h_phi, nper
      use vector_potentail_mod, only: ns, hs

      real(dp), parameter :: twopi = 2.d0*3.14159265358979d0

      real(dp), intent(in) :: s, theta, varphi
      real(dp), intent(out) :: ds, dtheta, dphi
      integer, intent(out) :: is, i_theta, i_phi

      ds = s/hs
      is = max(0, min(ns - 1, int(ds)))
      ds = (ds - dble(is))*hs
      is = is + 1

      dtheta = modulo(theta, twopi)/h_theta
      i_theta = max(0, min(n_theta - 1, int(dtheta)))
      dtheta = (dtheta - dble(i_theta))*h_theta
      i_theta = i_theta + 1

      dphi = modulo(varphi, twopi/dble(nper))/h_phi
      i_phi = max(0, min(n_phi - 1, int(dphi)))
      dphi = (dphi - dble(i_phi))*h_phi
      i_phi = i_phi + 1

   end subroutine normalize_coordinates

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine interpolate_vector_potential(s, ds, is, A_phi, A_theta, dA_phi_ds, dA_theta_ds)
      use vector_potentail_mod, only: torflux, sA_phi
      use new_vmec_stuff_mod, only: ns_A

      real(dp), intent(in) :: s, ds
      integer, intent(in) :: is
      real(dp), intent(out) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds
      integer :: k

      A_theta = torflux*s
      dA_theta_ds = torflux

      if (.not. allocated(sA_phi)) call spline_vmec_data
      A_phi = sA_phi(ns_A + 1, is)
      dA_phi_ds = 0.d0

      do k = ns_A, 1, -1
         A_phi = sA_phi(k, is) + ds*A_phi
         dA_phi_ds = sA_phi(k + 1, is)*dble(k) + ds*dA_phi_ds
      end do

   end subroutine interpolate_vector_potential

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine get_spline_coefficients_3d(is, i_theta, i_phi, stp_R, stp_Z, stp_lam, &
                                         dstp_R_ds, dstp_Z_ds, dstp_lam_ds)
      use new_vmec_stuff_mod, only: sR, sZ, slam, ns_s, ns_tp

      integer, intent(in) :: is, i_theta, i_phi
      real(dp), dimension(:, :), intent(out) :: stp_R, stp_Z, stp_lam
      real(dp), dimension(:, :), intent(out) :: dstp_R_ds, dstp_Z_ds, dstp_lam_ds
      integer :: nstp

      nstp = ns_tp + 1

      stp_R(1:nstp, 1:nstp) = sR(ns_s + 1, :, :, is, i_theta, i_phi)
      dstp_R_ds(1:nstp, 1:nstp) = 0.d0
      stp_Z(1:nstp, 1:nstp) = sZ(ns_s + 1, :, :, is, i_theta, i_phi)
      dstp_Z_ds(1:nstp, 1:nstp) = 0.d0
      stp_lam(1:nstp, 1:nstp) = slam(ns_s + 1, :, :, is, i_theta, i_phi)
      dstp_lam_ds(1:nstp, 1:nstp) = 0.d0

   end subroutine get_spline_coefficients_3d

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine interpolate_s_direction(ds, is, i_theta, i_phi, stp_R, stp_Z, stp_lam, &
                                      dstp_R_ds, dstp_Z_ds, dstp_lam_ds)
      use new_vmec_stuff_mod, only: sR, sZ, slam, ns_s, ns_tp

      real(dp), intent(in) :: ds
      integer, intent(in) :: is, i_theta, i_phi
      real(dp), dimension(:, :), intent(inout) :: stp_R, stp_Z, stp_lam
      real(dp), dimension(:, :), intent(inout) :: dstp_R_ds, dstp_Z_ds, dstp_lam_ds
      integer :: k, nstp, it, ip
      real(dp) :: ks

      nstp = ns_tp + 1

      do k = ns_s, 1, -1
         ks = real(k, dp)
         do ip = 1, nstp
            do it = 1, nstp
               stp_R(it, ip) = sR(k, it, ip, is, i_theta, i_phi) + &
                               ds*stp_R(it, ip)
               dstp_R_ds(it, ip) = sR(k + 1, it, ip, is, i_theta, i_phi)*ks + &
                                   ds*dstp_R_ds(it, ip)
               stp_Z(it, ip) = sZ(k, it, ip, is, i_theta, i_phi) + &
                               ds*stp_Z(it, ip)
               dstp_Z_ds(it, ip) = sZ(k + 1, it, ip, is, i_theta, i_phi)*ks + &
                                   ds*dstp_Z_ds(it, ip)
               stp_lam(it, ip) = slam(k, it, ip, is, i_theta, i_phi) + &
                                 ds*stp_lam(it, ip)
               dstp_lam_ds(it, ip) = slam(k + 1, it, ip, is, i_theta, i_phi)*ks + &
                                     ds*dstp_lam_ds(it, ip)
            end do
         end do
      end do

   end subroutine interpolate_s_direction

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine interpolate_theta_direction(dtheta, stp_R, stp_Z, stp_lam, dstp_R_ds, dstp_Z_ds, dstp_lam_ds, &
                                          sp_R, sp_Z, sp_lam, dsp_R_ds, dsp_Z_ds, dsp_lam_ds, &
                                          dsp_R_dt, dsp_Z_dt, dsp_lam_dt)
      use new_vmec_stuff_mod, only: ns_tp

      real(dp), intent(in) :: dtheta
      real(dp), dimension(:, :), intent(in) :: stp_R, stp_Z, stp_lam
      real(dp), dimension(:, :), intent(in) :: dstp_R_ds, dstp_Z_ds, dstp_lam_ds
      real(dp), dimension(:), intent(out) :: sp_R, sp_Z, sp_lam
      real(dp), dimension(:), intent(out) :: dsp_R_ds, dsp_Z_ds, dsp_lam_ds
      real(dp), dimension(:), intent(out) :: dsp_R_dt, dsp_Z_dt, dsp_lam_dt
      integer :: k, nstp

      nstp = ns_tp + 1

      sp_R(1:nstp) = stp_R(nstp, 1:nstp)
      dsp_R_ds(1:nstp) = dstp_R_ds(nstp, 1:nstp)
      dsp_R_dt(1:nstp) = 0.d0
      sp_Z(1:nstp) = stp_Z(nstp, 1:nstp)
      dsp_Z_ds(1:nstp) = dstp_Z_ds(nstp, 1:nstp)
      dsp_Z_dt(1:nstp) = 0.d0
      sp_lam(1:nstp) = stp_lam(nstp, 1:nstp)
      dsp_lam_ds(1:nstp) = dstp_lam_ds(nstp, 1:nstp)
      dsp_lam_dt(1:nstp) = 0.d0

      do k = ns_tp, 1, -1
         sp_R(1:nstp) = stp_R(k, 1:nstp) + dtheta*sp_R(1:nstp)
         dsp_R_ds(1:nstp) = dstp_R_ds(k, 1:nstp) + dtheta*dsp_R_ds(1:nstp)
         dsp_R_dt(1:nstp) = stp_R(k + 1, 1:nstp)*dble(k) + dtheta*dsp_R_dt(1:nstp)

         sp_Z(1:nstp) = stp_Z(k, 1:nstp) + dtheta*sp_Z(1:nstp)
         dsp_Z_ds(1:nstp) = dstp_Z_ds(k, 1:nstp) + dtheta*dsp_Z_ds(1:nstp)
         dsp_Z_dt(1:nstp) = stp_Z(k + 1, 1:nstp)*dble(k) + dtheta*dsp_Z_dt(1:nstp)

         sp_lam(1:nstp) = stp_lam(k, 1:nstp) + dtheta*sp_lam(1:nstp)
         dsp_lam_ds(1:nstp) = dstp_lam_ds(k, 1:nstp) + dtheta*dsp_lam_ds(1:nstp)
         dsp_lam_dt(1:nstp) = stp_lam(k + 1, 1:nstp)*dble(k) + dtheta*dsp_lam_dt(1:nstp)
      end do

   end subroutine interpolate_theta_direction

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine interpolate_phi_direction(dphi, sp_R, sp_Z, sp_lam, dsp_R_ds, dsp_Z_ds, dsp_lam_ds, &
                                        dsp_R_dt, dsp_Z_dt, dsp_lam_dt, &
                                        R, Z, alam, dR_ds, dZ_ds, dl_ds, dR_dt, dZ_dt, dl_dt, &
                                        dR_dp, dZ_dp, dl_dp)
      use new_vmec_stuff_mod, only: ns_tp

      real(dp), intent(in) :: dphi
      real(dp), dimension(:), intent(in) :: sp_R, sp_Z, sp_lam
      real(dp), dimension(:), intent(in) :: dsp_R_ds, dsp_Z_ds, dsp_lam_ds
      real(dp), dimension(:), intent(in) :: dsp_R_dt, dsp_Z_dt, dsp_lam_dt
      real(dp), intent(out) :: R, Z, alam, dR_ds, dZ_ds, dl_ds
      real(dp), intent(out) :: dR_dt, dZ_dt, dl_dt, dR_dp, dZ_dp, dl_dp
      integer :: k, nstp

      nstp = ns_tp + 1

      R = sp_R(nstp)
      dR_ds = dsp_R_ds(nstp)
      dR_dt = dsp_R_dt(nstp)
      dR_dp = 0.d0
      Z = sp_Z(nstp)
      dZ_ds = dsp_Z_ds(nstp)
      dZ_dt = dsp_Z_dt(nstp)
      dZ_dp = 0.d0
      alam = sp_lam(nstp)
      dl_ds = dsp_lam_ds(nstp)
      dl_dt = dsp_lam_dt(nstp)
      dl_dp = 0.d0

      do k = ns_tp, 1, -1
         R = sp_R(k) + dphi*R
         dR_ds = dsp_R_ds(k) + dphi*dR_ds
         dR_dt = dsp_R_dt(k) + dphi*dR_dt
         dR_dp = sp_R(k + 1)*dble(k) + dphi*dR_dp

         Z = sp_Z(k) + dphi*Z
         dZ_ds = dsp_Z_ds(k) + dphi*dZ_ds
         dZ_dt = dsp_Z_dt(k) + dphi*dZ_dt
         dZ_dp = sp_Z(k + 1)*dble(k) + dphi*dZ_dp

         alam = sp_lam(k) + dphi*alam
         dl_ds = dsp_lam_ds(k) + dphi*dl_ds
         dl_dt = dsp_lam_dt(k) + dphi*dl_dt
         dl_dp = sp_lam(k + 1)*dble(k) + dphi*dl_dp
      end do

   end subroutine interpolate_phi_direction

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine splint_vmec_data(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                               R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)
      use vector_potentail_mod, only: ns, hs

      real(dp), intent(in) :: s, theta, varphi
      real(dp), intent(out) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota
      real(dp), intent(out) :: R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
      real(dp), intent(out) :: dl_ds, dl_dt, dl_dp

      integer, parameter :: ns_max = 6
      integer :: is, i_theta, i_phi, is_rho
      real(dp) :: ds, dtheta, dphi, ds_rho, rho_tor
      real(dp), dimension(ns_max) :: sp_R, sp_Z, sp_lam
      real(dp), dimension(ns_max) :: dsp_R_ds, dsp_Z_ds, dsp_lam_ds
      real(dp), dimension(ns_max) :: dsp_R_dt, dsp_Z_dt, dsp_lam_dt
      real(dp), dimension(ns_max, ns_max) :: stp_R, stp_Z, stp_lam
      real(dp), dimension(ns_max, ns_max) :: dstp_R_ds, dstp_Z_ds, dstp_lam_ds

      call normalize_coordinates(s, theta, varphi, ds, is, dtheta, i_theta, dphi, i_phi)
      call interpolate_vector_potential(s, ds, is, A_phi, A_theta, dA_phi_ds, dA_theta_ds)

      aiota = -dA_phi_ds/dA_theta_ds
      rho_tor = sqrt(s)
      ds_rho = rho_tor/hs
      is_rho = max(0, min(ns - 2, int(ds_rho)))
      ds_rho = (ds_rho - dble(is_rho))*hs
      is_rho = is_rho + 1

      call get_spline_coefficients_3d(is_rho, i_theta, i_phi, stp_R, stp_Z, stp_lam, &
                                      dstp_R_ds, dstp_Z_ds, dstp_lam_ds)

      call interpolate_s_direction(ds_rho, is_rho, i_theta, i_phi, stp_R, stp_Z, stp_lam, &
                                   dstp_R_ds, dstp_Z_ds, dstp_lam_ds)

      call interpolate_theta_direction(dtheta, stp_R, stp_Z, stp_lam, &
                                       dstp_R_ds, dstp_Z_ds, dstp_lam_ds, &
                                       sp_R, sp_Z, sp_lam, dsp_R_ds, dsp_Z_ds, dsp_lam_ds, &
                                       dsp_R_dt, dsp_Z_dt, dsp_lam_dt)

      call interpolate_phi_direction(dphi, sp_R, sp_Z, sp_lam, &
                                     dsp_R_ds, dsp_Z_ds, dsp_lam_ds, &
                                     dsp_R_dt, dsp_Z_dt, dsp_lam_dt, &
                                     R, Z, alam, dR_ds, dZ_ds, dl_ds, &
                                     dR_dt, dZ_dt, dl_dt, dR_dp, dZ_dp, dl_dp)

      if (rho_tor > 0.0d0) then
         dR_ds = 0.5d0*dR_ds/rho_tor
         dZ_ds = 0.5d0*dZ_ds/rho_tor
         dl_ds = 0.5d0*dl_ds/rho_tor
      else
         ! Magnetic axis (s=0): d/ds carries the 1/(2*sqrt(s)) coordinate
         ! singularity and is undefined here. R, Z and lambda above are the exact
         ! axis position; callers that evaluate the axis use only those values.
         dR_ds = 0.0d0
         dZ_ds = 0.0d0
         dl_ds = 0.0d0
      end if

   end subroutine splint_vmec_data

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine splint_vmec_data_d2(s, theta, varphi, A_phi, A_theta, &
                                  dA_phi_ds, dA_theta_ds, aiota, R, Z, alam, &
                                  dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                  dl_ds, dl_dt, dl_dp, d2R, d2Z, d2l)
      !> Like splint_vmec_data but also returns the second derivatives of the
      !> coordinate map (R, Z, lambda) with respect to (s, vartheta, varphi).
      !> d2R, d2Z, d2l hold (ss, st, sp, tt, tp, pp) in that order.
      !>
      !> The 3D splines are stored in rho_tor = sqrt(s); first and second
      !> s-derivatives apply the chain rule for that substitution.
      use vector_potentail_mod, only: ns, hs
      use new_vmec_stuff_mod, only: sR, sZ, slam, ns_s, ns_tp

      real(dp), intent(in) :: s, theta, varphi
      real(dp), intent(out) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota
      real(dp), intent(out) :: R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
      real(dp), intent(out) :: dl_ds, dl_dt, dl_dp
      real(dp), intent(out) :: d2R(6), d2Z(6), d2l(6)

      integer :: is, i_theta, i_phi, is_rho
      integer :: ka, kb, kc, nrho, nt, np
      real(dp) :: ds, dtheta, dphi, ds_rho, rho_tor
      real(dp) :: pr, pt, pp, dpr, dpt, dpp, d2pr, d2pt, d2pp, w
      real(dp) :: fld(3), df(3, 3), d2f(3, 6)
      real(dp) :: drho, d2rho

      call normalize_coordinates(s, theta, varphi, ds, is, dtheta, i_theta, dphi, i_phi)
      call interpolate_vector_potential(s, ds, is, A_phi, A_theta, dA_phi_ds, dA_theta_ds)

      aiota = -dA_phi_ds/dA_theta_ds

      rho_tor = sqrt(s)
      ds_rho = rho_tor/hs
      is_rho = max(0, min(ns - 2, int(ds_rho)))
      ds_rho = (ds_rho - dble(is_rho))*hs
      is_rho = is_rho + 1

      nrho = ns_s + 1
      nt = ns_tp + 1
      np = ns_tp + 1

      fld = 0.d0
      df = 0.d0
      d2f = 0.d0

      do ka = 1, nrho
         pr = ds_rho**(ka - 1)
         dpr = 0.d0
         if (ka >= 2) dpr = dble(ka - 1)*ds_rho**(ka - 2)
         d2pr = 0.d0
         if (ka >= 3) d2pr = dble((ka - 1)*(ka - 2))*ds_rho**(ka - 3)

         do kb = 1, nt
            pt = dtheta**(kb - 1)
            dpt = 0.d0
            if (kb >= 2) dpt = dble(kb - 1)*dtheta**(kb - 2)
            d2pt = 0.d0
            if (kb >= 3) d2pt = dble((kb - 1)*(kb - 2))*dtheta**(kb - 3)

            do kc = 1, np
               pp = dphi**(kc - 1)
               dpp = 0.d0
               if (kc >= 2) dpp = dble(kc - 1)*dphi**(kc - 2)
               d2pp = 0.d0
               if (kc >= 3) d2pp = dble((kc - 1)*(kc - 2))*dphi**(kc - 3)

               call accumulate_map_term(sR(ka, kb, kc, is_rho, i_theta, i_phi), 1)
               call accumulate_map_term(sZ(ka, kb, kc, is_rho, i_theta, i_phi), 2)
               call accumulate_map_term(slam(ka, kb, kc, is_rho, i_theta, i_phi), 3)
            end do
         end do
      end do

      R = fld(1)
      Z = fld(2)
      alam = fld(3)

      ! Chain rule rho_tor = sqrt(s): d/ds = drho * d/drho with drho = 1/(2 rho).
      drho = 0.5d0/rho_tor
      d2rho = -0.25d0/(rho_tor**3)

      dR_dt = df(1, 2)
      dR_dp = df(1, 3)
      dZ_dt = df(2, 2)
      dZ_dp = df(2, 3)
      dl_dt = df(3, 2)
      dl_dp = df(3, 3)

      dR_ds = df(1, 1)*drho
      dZ_ds = df(2, 1)*drho
      dl_ds = df(3, 1)*drho

      ! d2 order (ss, st, sp, tt, tp, pp).
      d2R(1) = d2f(1, 1)*drho*drho + df(1, 1)*d2rho
      d2R(2) = d2f(1, 2)*drho
      d2R(3) = d2f(1, 3)*drho
      d2R(4) = d2f(1, 4)
      d2R(5) = d2f(1, 5)
      d2R(6) = d2f(1, 6)

      d2Z(1) = d2f(2, 1)*drho*drho + df(2, 1)*d2rho
      d2Z(2) = d2f(2, 2)*drho
      d2Z(3) = d2f(2, 3)*drho
      d2Z(4) = d2f(2, 4)
      d2Z(5) = d2f(2, 5)
      d2Z(6) = d2f(2, 6)

      d2l(1) = d2f(3, 1)*drho*drho + df(3, 1)*d2rho
      d2l(2) = d2f(3, 2)*drho
      d2l(3) = d2f(3, 3)*drho
      d2l(4) = d2f(3, 4)
      d2l(5) = d2f(3, 5)
      d2l(6) = d2f(3, 6)

   contains

      subroutine accumulate_map_term(coef, ifld)
         real(dp), intent(in) :: coef
         integer, intent(in) :: ifld

         w = coef
         fld(ifld) = fld(ifld) + w*pr*pt*pp

         df(ifld, 1) = df(ifld, 1) + w*dpr*pt*pp
         df(ifld, 2) = df(ifld, 2) + w*pr*dpt*pp
         df(ifld, 3) = df(ifld, 3) + w*pr*pt*dpp

         ! Second derivatives in rho/theta/phi, order (rr, rt, rp, tt, tp, pp).
         d2f(ifld, 1) = d2f(ifld, 1) + w*d2pr*pt*pp
         d2f(ifld, 2) = d2f(ifld, 2) + w*dpr*dpt*pp
         d2f(ifld, 3) = d2f(ifld, 3) + w*dpr*pt*dpp
         d2f(ifld, 4) = d2f(ifld, 4) + w*pr*d2pt*pp
         d2f(ifld, 5) = d2f(ifld, 5) + w*pr*dpt*dpp
         d2f(ifld, 6) = d2f(ifld, 6) + w*pr*pt*d2pp
      end subroutine accumulate_map_term

   end subroutine splint_vmec_data_d2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine splint_vmec_data_d3(s, theta, varphi, d3R, d3Z, d3l)
      !> Third derivatives of the coordinate map (R, Z, lambda) with respect to
      !> (s, vartheta, varphi). Packed order (10 components):
      !>   (sss, sst, ssp, stt, stp, spp, ttt, ttp, tpp, ppp).
      !> Needed for the analytic second derivative of the curvilinear metric
      !> (boozer_field_metric d2g_B pullback). Same rho_tor = sqrt(s) storage as
      !> splint_vmec_data_d2; the s-direction applies the rho->s chain rule to
      !> third order.
      use vector_potentail_mod, only: ns, hs
      use new_vmec_stuff_mod, only: sR, sZ, slam, ns_s, ns_tp

      real(dp), intent(in) :: s, theta, varphi
      real(dp), intent(out) :: d3R(10), d3Z(10), d3l(10)

      integer :: is, i_theta, i_phi, is_rho
      integer :: ka, kb, kc, nrho, nt, np
      real(dp) :: ds, dtheta, dphi, ds_rho, rho_tor
      real(dp) :: pr, pt, pp, dpr, dpt, dpp, d2pr, d2pt, d2pp, d3pr, d3pt, d3pp, w
      real(dp) :: df(3, 3), d2f(3, 6), d3f(3, 10)
      real(dp) :: drho, d2rho, d3rho

      call normalize_coordinates(s, theta, varphi, ds, is, dtheta, i_theta, dphi, i_phi)

      rho_tor = sqrt(s)
      ds_rho = rho_tor/hs
      is_rho = max(0, min(ns - 2, int(ds_rho)))
      ds_rho = (ds_rho - dble(is_rho))*hs
      is_rho = is_rho + 1

      nrho = ns_s + 1
      nt = ns_tp + 1
      np = ns_tp + 1

      df = 0.d0
      d2f = 0.d0
      d3f = 0.d0

      do ka = 1, nrho
         pr = ds_rho**(ka - 1)
         dpr = 0.d0
         if (ka >= 2) dpr = dble(ka - 1)*ds_rho**(ka - 2)
         d2pr = 0.d0
         if (ka >= 3) d2pr = dble((ka - 1)*(ka - 2))*ds_rho**(ka - 3)
         d3pr = 0.d0
         if (ka >= 4) d3pr = dble((ka - 1)*(ka - 2)*(ka - 3))*ds_rho**(ka - 4)

         do kb = 1, nt
            pt = dtheta**(kb - 1)
            dpt = 0.d0
            if (kb >= 2) dpt = dble(kb - 1)*dtheta**(kb - 2)
            d2pt = 0.d0
            if (kb >= 3) d2pt = dble((kb - 1)*(kb - 2))*dtheta**(kb - 3)
            d3pt = 0.d0
            if (kb >= 4) d3pt = dble((kb - 1)*(kb - 2)*(kb - 3))*dtheta**(kb - 4)

            do kc = 1, np
               pp = dphi**(kc - 1)
               dpp = 0.d0
               if (kc >= 2) dpp = dble(kc - 1)*dphi**(kc - 2)
               d2pp = 0.d0
               if (kc >= 3) d2pp = dble((kc - 1)*(kc - 2))*dphi**(kc - 3)
               d3pp = 0.d0
               if (kc >= 4) d3pp = dble((kc - 1)*(kc - 2)*(kc - 3))*dphi**(kc - 4)

               call accumulate_d3_term(sR(ka, kb, kc, is_rho, i_theta, i_phi), 1)
               call accumulate_d3_term(sZ(ka, kb, kc, is_rho, i_theta, i_phi), 2)
               call accumulate_d3_term(slam(ka, kb, kc, is_rho, i_theta, i_phi), 3)
            end do
         end do
      end do

      ! Chain rule rho_tor = sqrt(s): drho = 1/(2 rho), d2rho = -1/(4 rho^3),
      ! d3rho = 3/(8 rho^5).
      drho = 0.5d0/rho_tor
      d2rho = -0.25d0/(rho_tor**3)
      d3rho = 0.375d0/(rho_tor**5)

      call chain_d3(df(1, :), d2f(1, :), d3f(1, :), d3R)
      call chain_d3(df(2, :), d2f(2, :), d3f(2, :), d3Z)
      call chain_d3(df(3, :), d2f(3, :), d3f(3, :), d3l)

   contains

      subroutine accumulate_d3_term(coef, ifld)
         real(dp), intent(in) :: coef
         integer, intent(in) :: ifld

         w = coef
         ! rho first derivative and the rho-involving second derivatives that the
         ! s-direction chain rule references.
         df(ifld, 1) = df(ifld, 1) + w*dpr*pt*pp
         d2f(ifld, 1) = d2f(ifld, 1) + w*d2pr*pt*pp        ! rr
         d2f(ifld, 2) = d2f(ifld, 2) + w*dpr*dpt*pp        ! rt
         d2f(ifld, 3) = d2f(ifld, 3) + w*dpr*pt*dpp        ! rp

         ! Third derivatives in rho/theta/phi, order
         ! (rrr, rrt, rrp, rtt, rtp, rpp, ttt, ttp, tpp, ppp).
         d3f(ifld, 1) = d3f(ifld, 1) + w*d3pr*pt*pp
         d3f(ifld, 2) = d3f(ifld, 2) + w*d2pr*dpt*pp
         d3f(ifld, 3) = d3f(ifld, 3) + w*d2pr*pt*dpp
         d3f(ifld, 4) = d3f(ifld, 4) + w*dpr*d2pt*pp
         d3f(ifld, 5) = d3f(ifld, 5) + w*dpr*dpt*dpp
         d3f(ifld, 6) = d3f(ifld, 6) + w*dpr*pt*d2pp
         d3f(ifld, 7) = d3f(ifld, 7) + w*pr*d3pt*pp
         d3f(ifld, 8) = d3f(ifld, 8) + w*pr*d2pt*dpp
         d3f(ifld, 9) = d3f(ifld, 9) + w*pr*dpt*d2pp
         d3f(ifld, 10) = d3f(ifld, 10) + w*pr*pt*d3pp
      end subroutine accumulate_d3_term

      subroutine chain_d3(df1, d2f1, d3f1, out)
         !> Convert rho/theta/phi derivatives to s/theta/phi. df1 holds the rho
         !> first derivative in slot 1; d2f1 holds (rr, rt, rp) in slots 1..3;
         !> d3f1 holds the 10 rho-space third derivatives.
         real(dp), intent(in) :: df1(3), d2f1(6), d3f1(10)
         real(dp), intent(out) :: out(10)

         out(1) = d3f1(1)*drho**3 + 3.d0*drho*d2rho*d2f1(1) + d3rho*df1(1)  ! sss
         out(2) = d3f1(2)*drho*drho + d2rho*d2f1(2)                          ! sst
         out(3) = d3f1(3)*drho*drho + d2rho*d2f1(3)                          ! ssp
         out(4) = d3f1(4)*drho                                               ! stt
         out(5) = d3f1(5)*drho                                               ! stp
         out(6) = d3f1(6)*drho                                               ! spp
         out(7) = d3f1(7)                                                    ! ttt
         out(8) = d3f1(8)                                                    ! ttp
         out(9) = d3f1(9)                                                    ! tpp
         out(10) = d3f1(10)                                                  ! ppp
      end subroutine chain_d3

   end subroutine splint_vmec_data_d3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                         sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, &
                         Bcovar_s, Bcovar_vartheta, Bcovar_varphi)
      !> Evaluates A and B components in symmetry flux coordinates (s, vartheta, varphi)
      !> at a given point in VMEC coordinates (s, theta, varphi). theta is VMEC poloidal angle,
      !> vartheta is symmetry flux poloidal angle (also called theta* in VMEC)

      real(dp), intent(in) :: s, theta, varphi
      real(dp), intent(out) :: A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                               sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, &
                               Bcovar_s, Bcovar_vartheta, Bcovar_varphi
      real(dp) :: R, Z, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp

      call splint_vmec_data(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                            R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)

      call compute_field_components(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                    dA_theta_ds, dA_phi_ds, dl_ds, dl_dt, dl_dp, &
                                    sqg, Bctrvr_vartheta, Bctrvr_varphi, Bcovar_s, &
                                    Bcovar_vartheta, Bcovar_varphi)

   end subroutine vmec_field

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine compute_field_components(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                       dA_theta_ds, dA_phi_ds, dl_ds, dl_dt, dl_dp, &
                                       sqg, Bctrvr_vartheta, Bctrvr_varphi, Bcovar_s, &
                                       Bcovar_vartheta, Bcovar_varphi)

      real(dp), intent(in) :: R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                              dA_theta_ds, dA_phi_ds, dl_ds, dl_dt, dl_dp
      real(dp), intent(out) :: sqg, Bctrvr_vartheta, Bctrvr_varphi, Bcovar_vartheta, &
                               Bcovar_varphi, Bcovar_s

      real(dp) :: g(3, 3)

      call metric_tensor_symflux(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                 dl_ds, dl_dt, dl_dp, g, sqg)

      Bctrvr_vartheta = -dA_phi_ds/sqg
      Bctrvr_varphi = dA_theta_ds/sqg

      Bcovar_s = g(1, 2)*Bctrvr_vartheta + g(1, 3)*Bctrvr_varphi
      Bcovar_vartheta = g(2, 2)*Bctrvr_vartheta + g(2, 3)*Bctrvr_varphi
      Bcovar_varphi = g(3, 2)*Bctrvr_vartheta + g(3, 3)*Bctrvr_varphi
   end subroutine compute_field_components

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine metric_tensor_symflux(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                    dl_ds, dl_dt, dl_dp, g, sqg)
      !> Computes the metric tensor g and its square root sqg
      !> in symmetry flux coordinates (s, vartheta, varphi)
      !> at a given point in VMEC coordinates (s, theta, varphi).

      real(dp), intent(in) :: R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                              dl_ds, dl_dt, dl_dp
      real(dp), intent(out) :: g(3, 3), sqg

      real(dp), dimension(3, 3) :: cmat, gV
      real(dp) :: cjac, sqgV

      call metric_tensor_vmec(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, gV, sqgV)

      cjac = 1.d0/(1.d0 + dl_dt)
      sqg = sqgV*cjac

      cmat(1, 2:3) = 0.d0
      cmat(3, 1:2) = 0.d0
      cmat(1, 1) = 1.d0
      cmat(3, 3) = 1.d0
      cmat(2, 1) = -dl_ds*cjac
      cmat(2, 2) = cjac
      cmat(2, 3) = -dl_dp*cjac

      g = matmul(transpose(cmat), matmul(gV, cmat))
   end subroutine metric_tensor_symflux

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine christoffel_symflux(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
                                  dl_ds, dl_dt, dl_dp, d2R, d2Z, d2l, Gamma)
      !> Christoffel symbols of the second kind in symmetry flux coordinates
      !> (s, vartheta, varphi), Gamma(i,m,n) = Gamma^i_{mn}
      !>   = 0.5 g^{il} (d_m g_{ln} + d_n g_{lm} - d_l g_{mn}).
      !> The metric derivatives d_k g_{ij} are assembled algebraically from the
      !> first and second derivatives of R, Z and lambda, using the same
      !> g = C^T gV C construction as metric_tensor_symflux.

      real(dp), intent(in) :: R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
      real(dp), intent(in) :: dl_ds, dl_dt, dl_dp, d2R(6), d2Z(6), d2l(6)
      real(dp), intent(out) :: Gamma(3, 3, 3)

      integer :: i, m, n, l, k
      real(dp) :: g(3, 3), ginv(3, 3), dg(3, 3, 3), det
      real(dp) :: cmat(3, 3), gV(3, 3), dcmat(3, 3, 3), dgV(3, 3, 3)
      real(dp) :: cjac, dcjac(3), tmp(3, 3)
      ! First derivatives of R and Z grouped as dR(k) = dR/dx_k.
      real(dp) :: dR(3), dZ(3)
      ! Second derivatives as 3x3 symmetric matrices, d2R(i,j)=d^2 R/dx_i dx_j.
      real(dp) :: hR(3, 3), hZ(3, 3), hl(3, 3)
      integer :: idx6(3, 3)

      idx6 = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

      dR = [dR_ds, dR_dt, dR_dp]
      dZ = [dZ_ds, dZ_dt, dZ_dp]
      do i = 1, 3
         do k = 1, 3
            hR(i, k) = d2R(idx6(i, k))
            hZ(i, k) = d2Z(idx6(i, k))
            hl(i, k) = d2l(idx6(i, k))
         end do
      end do

      ! VMEC metric and its derivatives.
      gV(1, 1) = dR(1)**2 + dZ(1)**2
      gV(1, 2) = dR(1)*dR(2) + dZ(1)*dZ(2)
      gV(1, 3) = dR(1)*dR(3) + dZ(1)*dZ(3)
      gV(2, 2) = dR(2)**2 + dZ(2)**2
      gV(2, 3) = dR(2)*dR(3) + dZ(2)*dZ(3)
      gV(3, 3) = R**2 + dR(3)**2 + dZ(3)**2
      gV(2, 1) = gV(1, 2)
      gV(3, 1) = gV(1, 3)
      gV(3, 2) = gV(2, 3)

      do k = 1, 3
         dgV(1, 1, k) = 2.d0*(dR(1)*hR(1, k) + dZ(1)*hZ(1, k))
         dgV(1, 2, k) = hR(1, k)*dR(2) + dR(1)*hR(2, k) &
                        + hZ(1, k)*dZ(2) + dZ(1)*hZ(2, k)
         dgV(1, 3, k) = hR(1, k)*dR(3) + dR(1)*hR(3, k) &
                        + hZ(1, k)*dZ(3) + dZ(1)*hZ(3, k)
         dgV(2, 2, k) = 2.d0*(dR(2)*hR(2, k) + dZ(2)*hZ(2, k))
         dgV(2, 3, k) = hR(2, k)*dR(3) + dR(2)*hR(3, k) &
                        + hZ(2, k)*dZ(3) + dZ(2)*hZ(3, k)
         dgV(3, 3, k) = 2.d0*(dR(3)*hR(3, k) + dZ(3)*hZ(3, k))
         if (k == 1) dgV(3, 3, k) = dgV(3, 3, k) + 2.d0*R*dR(1)
         if (k == 2) dgV(3, 3, k) = dgV(3, 3, k) + 2.d0*R*dR(2)
         if (k == 3) dgV(3, 3, k) = dgV(3, 3, k) + 2.d0*R*dR(3)
         dgV(2, 1, k) = dgV(1, 2, k)
         dgV(3, 1, k) = dgV(1, 3, k)
         dgV(3, 2, k) = dgV(2, 3, k)
      end do

      ! Coordinate transform C from symflux to VMEC and its derivatives.
      cjac = 1.d0/(1.d0 + dl_dt)
      do k = 1, 3
         dcjac(k) = -cjac*cjac*hl(2, k)
      end do

      cmat = 0.d0
      cmat(1, 1) = 1.d0
      cmat(3, 3) = 1.d0
      cmat(2, 1) = -dl_ds*cjac
      cmat(2, 2) = cjac
      cmat(2, 3) = -dl_dp*cjac

      dcmat = 0.d0
      do k = 1, 3
         dcmat(2, 1, k) = -(hl(1, k)*cjac + dl_ds*dcjac(k))
         dcmat(2, 2, k) = dcjac(k)
         dcmat(2, 3, k) = -(hl(3, k)*cjac + dl_dp*dcjac(k))
      end do

      g = matmul(transpose(cmat), matmul(gV, cmat))

      ! d_k g = (d_k C)^T gV C + C^T (d_k gV) C + C^T gV (d_k C).
      do k = 1, 3
         tmp = matmul(transpose(dcmat(:, :, k)), matmul(gV, cmat)) &
               + matmul(transpose(cmat), matmul(dgV(:, :, k), cmat)) &
               + matmul(transpose(cmat), matmul(gV, dcmat(:, :, k)))
         dg(:, :, k) = tmp
      end do

      det = g(1, 1)*(g(2, 2)*g(3, 3) - g(2, 3)*g(3, 2)) &
            - g(1, 2)*(g(2, 1)*g(3, 3) - g(2, 3)*g(3, 1)) &
            + g(1, 3)*(g(2, 1)*g(3, 2) - g(2, 2)*g(3, 1))

      ginv(1, 1) = (g(2, 2)*g(3, 3) - g(2, 3)*g(3, 2))/det
      ginv(1, 2) = (g(1, 3)*g(3, 2) - g(1, 2)*g(3, 3))/det
      ginv(1, 3) = (g(1, 2)*g(2, 3) - g(1, 3)*g(2, 2))/det
      ginv(2, 1) = (g(2, 3)*g(3, 1) - g(2, 1)*g(3, 3))/det
      ginv(2, 2) = (g(1, 1)*g(3, 3) - g(1, 3)*g(3, 1))/det
      ginv(2, 3) = (g(1, 3)*g(2, 1) - g(1, 1)*g(2, 3))/det
      ginv(3, 1) = (g(2, 1)*g(3, 2) - g(2, 2)*g(3, 1))/det
      ginv(3, 2) = (g(1, 2)*g(3, 1) - g(1, 1)*g(3, 2))/det
      ginv(3, 3) = (g(1, 1)*g(2, 2) - g(1, 2)*g(2, 1))/det

      Gamma = 0.d0
      do i = 1, 3
         do m = 1, 3
            do n = 1, 3
               do l = 1, 3
                  Gamma(i, m, n) = Gamma(i, m, n) + 0.5d0*ginv(i, l)* &
                                   (dg(l, n, m) + dg(l, m, n) - dg(m, n, l))
               end do
            end do
         end do
      end do

   end subroutine christoffel_symflux

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine metric_tensor_vmec(R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, gV, sqgV)
      !> Computes the metric tensor g and its square root sqgV
      !> in VMEC coordinates (s, theta, varphi).

      real(dp), intent(in) :: R, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp
      real(dp), intent(out) :: gV(3, 3), sqgV

      gV(1, 1) = dR_ds**2 + dZ_ds**2
      gV(1, 2) = dR_ds*dR_dt + dZ_ds*dZ_dt
      gV(1, 3) = dR_ds*dR_dp + dZ_ds*dZ_dp
      gV(2, 1) = gV(1, 2)
      gV(2, 2) = dR_dt**2 + dZ_dt**2
      gV(2, 3) = dR_dt*dR_dp + dZ_dt*dZ_dp
      gV(3, 1) = gV(1, 3)
      gV(3, 2) = gV(2, 3)
      gV(3, 3) = R**2 + dR_dp**2 + dZ_dp**2
      sqgV = R*(dR_dt*dZ_ds - dR_ds*dZ_dt)
   end subroutine metric_tensor_vmec

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine splint_iota(s, aiota, daiota_ds)

      use vector_potentail_mod, only: ns, hs, torflux, sA_phi
      use new_vmec_stuff_mod, only: ns_A

      integer :: is, k
      real(dp) :: s, ds, dA_phi_ds, dA_theta_ds, d2A_phi_ds2, aiota, daiota_ds

      dA_theta_ds = torflux

      ds = s/hs
      is = max(0, min(ns - 1, int(ds)))
      ds = (ds - dble(is))*hs
      is = is + 1

      dA_phi_ds = 0.d0

      do k = ns_A, 1, -1
         dA_phi_ds = sA_phi(k + 1, is)*dble(k) + ds*dA_phi_ds
      end do

      d2A_phi_ds2 = 0.d0

      do k = ns_A, 2, -1
         d2A_phi_ds2 = sA_phi(k + 1, is)*dble(k)*dble(k - 1) + ds*d2A_phi_ds2
      end do

      aiota = -dA_phi_ds/dA_theta_ds
      daiota_ds = -d2A_phi_ds2/dA_theta_ds

   end subroutine splint_iota

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine splint_lambda(s, theta, varphi, alam, dl_dt)

      use new_vmec_stuff_mod, only: n_theta, n_phi, h_theta, h_phi, slam, nper, ns_s, ns_tp
      use vector_potentail_mod, only: ns, hs

      real(dp), parameter :: twopi = 2.d0*3.14159265358979d0

      integer :: is, i_theta, i_phi, k
      real(dp) :: ds, dtheta, dphi
      real(dp) :: s, theta, varphi, alam, dl_dt

      integer, parameter :: ns_max = 6

      integer :: nstp

      real(dp), dimension(ns_max) :: sp_lam
      real(dp), dimension(ns_max) :: dsp_lam_dt
      real(dp), dimension(ns_max, ns_max) :: stp_lam

      nstp = ns_tp + 1

      ds = sqrt(s)/hs
      is = max(0, min(ns - 1, int(ds)))
      ds = (ds - dble(is))*hs
      is = is + 1

      dtheta = modulo(theta, twopi)/h_theta
      i_theta = max(0, min(n_theta - 1, int(dtheta)))
      dtheta = (dtheta - dble(i_theta))*h_theta
      i_theta = i_theta + 1

      dphi = modulo(varphi, twopi/dble(nper))/h_phi
      i_phi = max(0, min(n_phi - 1, int(dphi)))
      dphi = (dphi - dble(i_phi))*h_phi
      i_phi = i_phi + 1

      ! Begin interpolation over $s$

      stp_lam(1:nstp, 1:nstp) = slam(ns_s + 1, :, :, is, i_theta, i_phi)

      do k = ns_s, 1, -1
         stp_lam(1:nstp, 1:nstp) = slam(k, :, :, is, i_theta, i_phi) + ds*stp_lam(1:nstp, 1:nstp)
      end do

      ! End interpolation over $s$
      !----------------------------

      ! Begin interpolation over $\theta$

      sp_lam(1:nstp) = stp_lam(nstp, 1:nstp)
      dsp_lam_dt(1:nstp) = 0.d0

      do k = ns_tp, 1, -1
         sp_lam(1:nstp) = stp_lam(k, 1:nstp) + dtheta*sp_lam(1:nstp)
         dsp_lam_dt(1:nstp) = stp_lam(k + 1, 1:nstp)*dble(k) + dtheta*dsp_lam_dt(1:nstp)
      end do

      ! End interpolation over $\theta$
      !--------------------------------

      ! Begin interpolation over $\varphi$

      alam = sp_lam(nstp)
      dl_dt = dsp_lam_dt(nstp)

      do k = ns_tp, 1, -1
         alam = sp_lam(k) + dphi*alam
         dl_dt = dsp_lam_dt(k) + dphi*dl_dt
      end do

   end subroutine splint_lambda

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   ! Go from s to rho grid, with special treatment of the axis.

   ! Interpolate values from s to rho grid. It is assumed that the
   ! innermost points of the input grid are not valid/noisy and thus need
   ! special treatment.
   ! This is done by extrapolating from values outside of this region to
   ! the axis.
   ! Extrapolation is done linear (more robust).
   ! An intermediate rescaling with rho can be used (might be useful to
   ! enforce behaviour near the axis). This will be in effect for
   ! extrapolating to the axis and for the interpolation to the new grid.

   ! input:
   ! ------
   ! m: integer, exponent, values <= 0 are ignored. Intermediate scaling of
   !   values is done with rho**m.
   ! ns: integer, size of input array.
   ! nrho: integer, size of output array.
   ! nheal: integer,
   ! arr_in: real(dp) 1d array, with ns elements.

   ! output:
   ! -------
   ! arr_out: real(dp) 1d array, with nrho elements.

   ! sideeffects:
   ! ------------
   ! none
   subroutine s_to_rho_healaxis(m, ns, nrho, nheal, arr_in, arr_out)

      use new_vmec_stuff_mod, only: ns_s, old_axis_healing

      integer, intent(in) :: m, ns, nrho, nheal
      real(dp), dimension(ns), intent(in) :: arr_in
      real(dp), dimension(nrho), intent(out) :: arr_out

      integer :: irho, is, k, nhe
      real(dp) :: hs, hrho, s, ds, rho, a, b, c
      real(dp), dimension(:, :), allocatable :: splcoe

      hs = 1.d0/dble(ns - 1)
      hrho = 1.d0/dble(nrho - 1)

      nhe = max(1, nheal) + 1

      ! Rescale
      do is = nhe, ns
         if (m .gt. 0) then
            rho = sqrt(hs*dble(is - 1))
            arr_out(is) = arr_in(is)/rho**m
         else
            arr_out(is) = arr_in(is)
         end if
      end do

      if (old_axis_healing) then
         ! parabolic extrapolation:
         a = arr_out(nhe)
         b = 0.5d0*(4.d0*arr_out(nhe + 1) - 3.d0*arr_out(nhe) - arr_out(nhe + 2))
         c = 0.5d0*(arr_out(nhe) + arr_out(nhe + 2) - 2.d0*arr_out(nhe + 1))

         do is = 1, nhe - 1
            arr_out(is) = a + b*dble(is - nhe) + c*dble(is - nhe)**2
         end do

      else
         ! linear extrapolation ("less accurate" but more robust):
         a = arr_out(nhe)
         b = arr_out(nhe + 1) - arr_out(nhe)

         do is = 1, nhe - 1
            arr_out(is) = a + b*dble(is - nhe)
         end do

      end if

      allocate (splcoe(0:ns_s, ns))

      splcoe(0, :) = arr_out

      call spl_reg(ns_s, ns, hs, splcoe)

      do irho = 1, nrho
         rho = hrho*dble(irho - 1)
         s = rho**2

         ds = s/hs
         is = max(0, min(ns - 1, int(ds)))
         ds = (ds - dble(is))*hs
         is = is + 1

         arr_out(irho) = splcoe(ns_s, is)

         do k = ns_s - 1, 0, -1
            arr_out(irho) = splcoe(k, is) + ds*arr_out(irho)
         end do

         ! Undo rescaling
         if (m .gt. 0) arr_out(irho) = arr_out(irho)*rho**m
      end do

      deallocate (splcoe)

   end subroutine s_to_rho_healaxis

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine s_to_rho_polyfit(m, ns, nrho, i_anchor, ndeg, arr_in, arr_out, anchor)
      !> Resamples one Fourier amplitude from the uniform s grid to the uniform
      !> rho = sqrt(s) grid, enforcing analytic axis regularity without the noise
      !> amplification of the pure power-law continuation. The rescaled amplitude
      !> g(s) = c(s)/rho**|m| is even-analytic in rho, hence smooth in s. Surfaces
      !> inside i_anchor are unreliable; g there is replaced by a least-squares
      !> polynomial in s (a Zernike radial polynomial once rho**|m| is restored)
      !> fitted to a window of reliable surfaces just outside i_anchor. g is then
      !> splined over the full grid and rho**|m| restored: c(rho) = rho**|m|*P(s),
      !> exact rho**|m| at the axis, smooth across the match, reproducing the
      !> reliable surfaces.
      use new_vmec_stuff_mod, only: ns_s

      integer, intent(in) :: m, ns, nrho, i_anchor, ndeg
      real(dp), dimension(ns), intent(in) :: arr_in
      real(dp), dimension(nrho), intent(out) :: arr_out
      !> anchor: exact s=0 value of the m=0 amplitude (the magnetic axis,
      !> raxis_cc/zaxis_cs). When present and m=0, the inward extrapolation is
      !> pinned to it so the healed interior meets the exact axis with no kink.
      real(dp), intent(in), optional :: anchor

      integer, parameter :: m_clamp = 50
      integer :: irho, is, k, mc, nwin, d, j
      real(dp) :: hs, hrho, s, ds, rho, u, s_c, s_scale, s_anchor, delta
      real(dp), dimension(:), allocatable :: g, swin, gwin, coef
      real(dp), dimension(:, :), allocatable :: splcoe

      hs = 1.d0/dble(ns - 1)
      hrho = 1.d0/dble(nrho - 1)
      mc = min(m, m_clamp)

      allocate (g(ns))
      if (mc > 0) then
         do is = i_anchor, ns
            rho = sqrt(hs*dble(is - 1))
            g(is) = arr_in(is)/rho**mc
         end do
      else
         g(i_anchor:ns) = arr_in(i_anchor:ns)
      end if

      ! Least-squares polynomial in s through a window of reliable surfaces,
      ! used to extrapolate g smoothly inward to the axis.
      d = max(0, min(ndeg, ns - i_anchor))
      nwin = min(ns - i_anchor + 1, 2*(d + 1) + 2)
      allocate (swin(nwin), gwin(nwin), coef(0:d))
      do j = 1, nwin
         swin(j) = hs*dble(i_anchor - 1 + j - 1)
         gwin(j) = g(i_anchor + j - 1)
      end do
      call polyfit_s(swin, gwin, nwin, d, coef, s_c, s_scale)
      do is = 1, i_anchor - 1
         s = hs*dble(is - 1)
         u = (s - s_c)/s_scale
         g(is) = coef(d)
         do k = d - 1, 0, -1
            g(is) = coef(k) + u*g(is)
         end do
      end do

      ! Pin the m=0 axis amplitude to the exact magnetic axis value, blended
      ! smoothly (vanishes at the reliable-window edge i_anchor) so the healed
      ! interior reaches the exact axis at s=0 with no kink.
      if (mc == 0 .and. present(anchor) .and. i_anchor > 1) then
         s_anchor = hs*dble(i_anchor - 1)
         delta = anchor - g(1)
         do is = 1, i_anchor - 1
            s = hs*dble(is - 1)
            g(is) = g(is) + delta*(1.0_dp - s/s_anchor)**2
         end do
      end if

      allocate (splcoe(0:ns_s, ns))
      splcoe(0, :) = g
      call spl_reg(ns_s, ns, hs, splcoe)
      do irho = 1, nrho
         rho = hrho*dble(irho - 1)
         s = rho**2
         ds = s/hs
         is = max(0, min(ns - 1, int(ds)))
         ds = (ds - dble(is))*hs
         is = is + 1
         arr_out(irho) = splcoe(ns_s, is)
         do k = ns_s - 1, 0, -1
            arr_out(irho) = splcoe(k, is) + ds*arr_out(irho)
         end do
         if (mc > 0) arr_out(irho) = arr_out(irho)*rho**mc
      end do

      deallocate (g, swin, gwin, coef, splcoe)
   end subroutine s_to_rho_polyfit

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine polyfit_s(s, g, n, d, coef, s_c, s_scale)
      !> Least-squares fit g(j) ~ sum_k coef(k)*u**k with u=(s-s_c)/s_scale,
      !> degree d, via the normal equations. Centering and scaling keep the
      !> Vandermonde well conditioned for the small near-axis s window.
      integer, intent(in) :: n, d
      real(dp), intent(in) :: s(n), g(n)
      real(dp), intent(out) :: coef(0:d), s_c, s_scale

      integer :: i, k, l
      real(dp) :: u
      real(dp), allocatable :: mat(:, :), rhs(:), upow(:)

      s_c = sum(s)/dble(n)
      s_scale = 0.5d0*(maxval(s) - minval(s))
      if (s_scale <= 0.d0) s_scale = 1.d0

      allocate (mat(0:d, 0:d), rhs(0:d), upow(0:2*d))
      mat = 0.d0
      rhs = 0.d0
      do i = 1, n
         u = (s(i) - s_c)/s_scale
         upow(0) = 1.d0
         do k = 1, 2*d
            upow(k) = upow(k - 1)*u
         end do
         do k = 0, d
            do l = 0, d
               mat(k, l) = mat(k, l) + upow(k + l)
            end do
            rhs(k) = rhs(k) + upow(k)*g(i)
         end do
      end do
      call solve_small(mat, rhs, coef, d + 1)
      deallocate (mat, rhs, upow)
   end subroutine polyfit_s

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine solve_small(a, b, x, n)
      !> Dense solve a x = b (n small) by Gaussian elimination with partial
      !> pivoting. a and b are consumed.
      integer, intent(in) :: n
      real(dp), intent(inout) :: a(n, n), b(n)
      real(dp), intent(out) :: x(n)

      integer :: i, j, k, ip
      real(dp) :: piv, factor, tmp

      do k = 1, n - 1
         ip = k
         piv = abs(a(k, k))
         do i = k + 1, n
            if (abs(a(i, k)) > piv) then
               piv = abs(a(i, k))
               ip = i
            end if
         end do
         if (ip /= k) then
            do j = k, n
               tmp = a(k, j); a(k, j) = a(ip, j); a(ip, j) = tmp
            end do
            tmp = b(k); b(k) = b(ip); b(ip) = tmp
         end if
         do i = k + 1, n
            factor = a(i, k)/a(k, k)
            do j = k, n
               a(i, j) = a(i, j) - factor*a(k, j)
            end do
            b(i) = b(i) - factor*b(k)
         end do
      end do
      do i = n, 1, -1
         x(i) = b(i)
         do j = i + 1, n
            x(i) = x(i) - a(i, j)*x(j)
         end do
         x(i) = x(i)/a(i, i)
      end do
   end subroutine solve_small

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   function axis_anchor_index(ns) result(i_anchor)
      !> Innermost reliable full-grid surface: the first surface at or
      !> outside s_axis_heal. VMEC full-grid harmonics inside this
      !> startup layer can violate the rho**|m| analyticity condition, so
      !> they are discarded and replaced by the selected continuation.
      use new_vmec_stuff_mod, only: s_axis_heal, rho_axis_heal, ns_s, &
                                    rho_axis_heal_warning_printed

      integer, intent(in) :: ns
      integer :: i_anchor
      real(dp) :: s_heal

      if (rho_axis_heal > 0.0d0) then
         s_heal = rho_axis_heal**2
         if (.not. rho_axis_heal_warning_printed) then
            print *, 'WARNING: rho_axis_heal is deprecated; use s_axis_heal = ', s_heal
            print *, '         This setting will become an error in the next SIMPLE release.'
            rho_axis_heal_warning_printed = .True.
         end if
      else
         s_heal = s_axis_heal
      end if

      i_anchor = 1 + nint(s_heal*dble(ns - 1))
      i_anchor = max(2, min(i_anchor, ns - ns_s - 1))
   end function axis_anchor_index

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine s_to_rho_power_law(m, ns, nrho, i_anchor, arr_in, arr_out)
      !> Resamples one Fourier amplitude from the uniform s grid to the
      !> uniform rho = sqrt(s) grid, enforcing analytic regularity at the
      !> axis. Surfaces below i_anchor are discarded: the amplitude is
      !> continued to the axis as c(rho) = c_anchor*(rho/rho_anchor)**|m|
      !> (m clamped at 50), anchored at the innermost reliable surface.
      !> Outside the anchor the spline acts on c/rho**|m|, an even analytic
      !> function of rho, built only from the reliable surfaces. This is the
      !> booz_xform_to_boozer_chartmap continuation applied at the VMEC
      !> harmonic level, so both the Boozer and canonical transforms inherit
      !> a clean near-axis field from one place.
      use new_vmec_stuff_mod, only: ns_s

      integer, intent(in) :: m, ns, nrho, i_anchor
      real(dp), dimension(ns), intent(in) :: arr_in
      real(dp), dimension(nrho), intent(out) :: arr_out

      integer, parameter :: m_clamp = 50

      integer :: irho, is, k, mc, nsub
      real(dp) :: hs, hrho, s, s_anchor, ds, rho, rho_anchor
      real(dp), dimension(:, :), allocatable :: splcoe

      hs = 1.d0/dble(ns - 1)
      hrho = 1.d0/dble(nrho - 1)
      mc = min(m, m_clamp)
      s_anchor = hs*dble(i_anchor - 1)
      rho_anchor = sqrt(s_anchor)
      nsub = ns - i_anchor + 1

      allocate (splcoe(0:ns_s, nsub))

      if (mc > 0) then
         do is = i_anchor, ns
            rho = sqrt(hs*dble(is - 1))
            splcoe(0, is - i_anchor + 1) = arr_in(is)/rho**mc
         end do
      else
         splcoe(0, :) = arr_in(i_anchor:ns)
      end if

      call spl_reg(ns_s, nsub, hs, splcoe)

      do irho = 1, nrho
         rho = hrho*dble(irho - 1)
         if (rho < rho_anchor) then
            arr_out(irho) = splcoe(0, 1)*rho**mc
         else
            s = rho**2
            ds = (s - s_anchor)/hs
            is = max(0, min(nsub - 1, int(ds)))
            ds = (ds - dble(is))*hs
            is = is + 1
            arr_out(irho) = splcoe(ns_s, is)
            do k = ns_s - 1, 0, -1
               arr_out(irho) = splcoe(k, is) + ds*arr_out(irho)
            end do
            if (mc > 0) arr_out(irho) = arr_out(irho)*rho**mc
         end if
      end do

      deallocate (splcoe)
   end subroutine s_to_rho_power_law

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine determine_nheal_for_axis(m, ns, arr_in, nheal)
      !> Determines the number of first radial points, nheal, where data is
      !> replaced by extrapolation.
      !>
      !> Takes cosine harmonic amplitude of arr_in and checks the difference
      !> between the data value at point i and the value obtained by 3-rd
      !> order Lagrange polynomial extrapolation from the points i+1,i+2,i+3
      !> and i+4. If relative value exceeds given tolerance (parameter set
      !> to 30%) this point is regarded as "bad" so that amplitude values
      !> for all points with indices smaller than i are replaced by linear
      !> extrapolation from points i and i+1. Note that harmonic function at
      !> 60 points per period is extrapolated with relative error about
      !> 1e-4, therefore 30% is quite a large tolerance which is exceeded if
      !> there is noise in the data. Even with this high tolerance, some
      !> harmonics have to be extrapolated almost from the very edge. Those,
      !> however, have amplitudes about 8 orders of magnitude smaller than
      !> main harmonics and, therefore, play no role.
      !>
      !> input:
      !> ------
      !> m:
      !> ns:
      !> arr_in: real(dp) array (ns entries), data from which to
      !>   determine the number of points to extrapolate at the axis.
      !>
      !> output:
      !> -------
      !> nheal: integer, number of points to extrapolate at the axis.

      ! Lagrange polynomial stencil size for checking the data by extraplation:
      integer, parameter :: nplag = 4
      ! tolerance for Lagrange polynomial extrapolation by one point (to check if data is noisy):
      real(dp), parameter :: tol = 3.d-1
      real(dp), parameter :: tiny = 1.d-200
      ! 3-rd order Lagrange polynomial extrapolation coefficients from points (1,2,3,4) to point 0:
      real(dp), parameter, dimension(nplag) :: weight = (/4.d0, -6.d0, 4.d0, -1.d0/)

      integer, intent(in) :: m, ns
      integer, intent(out) :: nheal
      real(dp), dimension(ns), intent(in) :: arr_in

      integer :: is, ncheck

      real(dp) :: hs, rho, rho_nonzero, errmax

      real(dp), dimension(:), allocatable :: arr

      ! We check points which are away by more than 3 stencils from the edge:
      ncheck = ns - 3*nplag

      hs = 1.d0/dble(ns - 1)
      allocate (arr(ns))

      do is = 2, ns
         if (m > 0) then
            rho = sqrt(hs*dble(is - 1))
            rho_nonzero = max(rho**m, tiny)
            arr(is) = arr_in(is)/rho_nonzero
         else
            arr(is) = arr_in(is)
         end if
      end do

      nheal = 1
      do is = ncheck, 2, -1
         nheal = is
         errmax = maxval(abs(arr(is:is + nplag)))*tol
         if (abs(arr(is) - sum(arr(is + 1:is + nplag)*weight)) > errmax) then
            exit
         end if
      end do

      deallocate (arr)

   end subroutine determine_nheal_for_axis

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine volume_and_B00(volume, B00)

      use new_vmec_stuff_mod, only: n_theta, n_phi, h_theta, h_phi, nper

      integer :: i_theta, i_phi
      real(dp) :: volume, B00
      real(dp) :: B3, B2, bmod2
      real(dp) :: s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                  R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp, &
                  sqg, Bctrvr_vartheta, Bctrvr_varphi, &
                  Bcovar_s, Bcovar_vartheta, Bcovar_varphi

      s = 0.9999999999d0
      volume = 0.d0

      do i_theta = 0, n_theta - 2
         theta = h_theta*dble(i_theta)
         do i_phi = 0, n_phi - 2
            varphi = h_phi*dble(i_phi)

            call splint_vmec_data(s, theta, varphi, A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, &
                                  R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp)

            volume = volume + R**2*dZ_dt
         end do
      end do

      volume = 0.5d0*abs(volume)*h_theta*h_phi*dble(nper)

      s = 1d-8
      theta = 0.d0
      B2 = 0.d0
      B3 = 0.d0
      do i_phi = 0, n_phi - 2
         varphi = h_phi*dble(i_phi)

         call vmec_field(s, theta, varphi, A_theta, A_phi, dA_theta_ds, dA_phi_ds, aiota, &
                         sqg, alam, dl_ds, dl_dt, dl_dp, Bctrvr_vartheta, Bctrvr_varphi, &
                         Bcovar_s, Bcovar_vartheta, Bcovar_varphi)

         bmod2 = Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi
         B2 = B2 + bmod2/Bctrvr_varphi
         B3 = B3 + bmod2*sqrt(bmod2)/Bctrvr_varphi
      end do

      B00 = B3/B2

   end subroutine volume_and_B00
end module spline_vmec_sub
