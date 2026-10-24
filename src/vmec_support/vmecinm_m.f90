module vmecin_sub
    implicit none
contains
    subroutine vmecin(rmnc, zmns, almns, rmns, zmnc, almnc, aiota, phi, sps, axm, axn, s, &
                      nsurfb, nstrb, kparb, flux)
        !-----
        ! Usage:
        !
        !    call vmecin(rmnc,zmns,almns,rmns,zmnc,almnc,aiota,phi,sps,axm,axn,s,    &
        !               nsurfm,nstrm,kpar,torflux)
        !
        !  where scalars are:
        !
        !  nstrm - (integer) number of harmonics
        !  kpar  - (integer) number of radial points
        !
        !  vectors are:
        !  aiota(0:kpar)             - (double precision) iota profile
        !  sps(0:kpar) = 0:kpar      - (double precision) radial index as dble number
        !  phi(0:kpar)               - (double precision) toroidal flux
        !  s(0:kpar)                 - (double precision) normalized toroidal flux
        !  axm(nstrm)                - (double precision) poloidal mode numbers
        !  axn(nstrm)                - (double precision) toroidal mode numbers
        !
        !  matrices are:
        !  rmnc(nstrm,0:kpar)     - (double precision) profiles of Fourier
        !  amplitudes of R for various cos(m*theta - n*phi) harmonics
        !  zmnc(nstrm,0:kpar)     - (double precision) profiles of Fourier
        !  amplitudes of Z for various cos(m*theta - n*phi) harmonics
        !  almnc(nstrm,0:kpar)     - (double precision) profiles of Fourier
        !  amplitudes of lambda for various cos(m*theta - n*phi) harmonics
        !  rmns(nstrm,0:kpar)     - (double precision) profiles of Fourier
        !  amplitudes of R for various sin(m*theta - n*phi) harmonics
        !  zmns(nstrm,0:kpar)     - (double precision) profiles of Fourier
        !  amplitudes of Z for various sin(m*theta - n*phi) harmonics
        !  almns(nstrm,0:kpar)     - (double precision) profiles of Fourier
        !  amplitudes of lambda for various sin(m*theta - n*phi) harmonics

        use libneo_kinds, only: dp
        !use math_constants, only: PI
        use new_vmec_stuff_mod, only: netcdffile, vmec_B_scale, vmec_RZ_scale, rmajor, &
                                      raxis_cc, zaxis_cs, raxis_cs, zaxis_cc
        use nctools_module, only: nc_open, nc_close, nc_get

        real(dp), parameter :: fac_b0 = 1d4, fac_r0 = 1d2, EPS = 1d-10, pi=3.14159265358979d0
        ! TODO: Replace pi by math_constants more accurate version. This requires new golden record in SIMPLE
        real(dp) :: fac_b, fac_r

        integer :: nsurfb, nstrb, kparb, ncid, i
        real(dp) :: flux
        real(dp), dimension(nstrb) :: axm, axn
        real(dp), dimension(0:kparb) :: sps, aiota, phi, s
        real(dp), dimension(nstrb, 0:kparb) :: rmnc, zmnc, almnc, lmnc
        real(dp), dimension(nstrb, 0:kparb) :: rmns, zmns, almns, lmns
        integer :: lasym_int
        logical :: lasym
        character(len=1000) :: filename

        associate (dummy => nsurfb)
        end associate

        fac_b = fac_b0*vmec_B_scale
        fac_r = fac_r0*vmec_RZ_scale

        do i = 0, kparb
            sps(i) = dble(i)
        end do

        filename = netcdffile
        if (len_trim(filename) == 0) filename = 'wout.nc'
        call nc_open(trim(filename), ncid)

        call nc_get(ncid, 'lasym__logical__', lasym_int)
        lasym = (lasym_int == 1)

        call nc_get(ncid, 'phi', phi)
        ! VMEC zeta is the physical geometric cylindrical azimuth.  Its usual
        ! (s, theta, zeta) chart is nevertheless left-handed (signgs=-1), because
        ! of the poloidal-angle orientation.  For an oriented Jacobian J<0,
        ! B^zeta=(dA_theta/ds)/J; hence a positive outward VMEC Phi_tor requires
        ! A_theta=-Phi_tor/(2*pi).  This is why SIMPLE stores
        ! torflux=-phi_vmec/(2*pi): no toroidal-angle reversal is involved.
        ! A_phi=-chi then gives iota=-dA_phi/dA_theta=+iota_vmec.
        phi = -phi/(2*pi)

        flux = phi(kparb)
        flux = flux*fac_b*fac_r**2
        phi = phi*fac_b*fac_r**2
        s = phi/flux

        call nc_get(ncid, 'xm', axm)
        call nc_get(ncid, 'xn', axn)
        call nc_get(ncid, 'iotaf', aiota)
        call nc_get(ncid, 'rmnc', rmnc)
        call nc_get(ncid, 'zmns', zmns)
        call nc_get(ncid, 'lmns', lmns)
        ! Exact magnetic axis harmonics (n=0..ntor), independent of the flux
        ! surfaces. Used for the analytic near-axis chartmap limit at s=0.
        if (allocated(raxis_cc)) then
            call nc_get(ncid, 'raxis_cc', raxis_cc)
            call nc_get(ncid, 'zaxis_cs', zaxis_cs)
        end if
        if (lasym) then
            call nc_get(ncid, 'rmns', rmns)
            call nc_get(ncid, 'zmnc', zmnc)
            call nc_get(ncid, 'lmnc', lmnc)
            if (allocated(raxis_cs)) then
                call nc_get(ncid, 'raxis_cs', raxis_cs)
                call nc_get(ncid, 'zaxis_cc', zaxis_cc)
            end if
        else
            rmns = 0d0
            zmnc = 0d0
            lmnc = 0d0
            if (allocated(raxis_cs)) then
                raxis_cs = 0d0
                zaxis_cc = 0d0
            end if
        end if
        ! Convert half-mesh to full mesh for lambda
        ! added by Christopher Albert, 2020-02-11
        call lambda_half_to_full_map(lmnc, almnc, kparb)
        call lambda_half_to_full_map(lmns, almns, kparb)

! Fallback if rmajor defined by volume is not given
        if (abs(rmajor - 0d0) < EPS) rmajor = rmnc(1, 0)
        rmnc = rmnc*fac_r
        zmnc = zmnc*fac_r
        rmns = rmns*fac_r
        zmns = zmns*fac_r
        if (allocated(raxis_cc)) then
            raxis_cc = raxis_cc*fac_r
            zaxis_cs = zaxis_cs*fac_r
            raxis_cs = raxis_cs*fac_r
            zaxis_cc = zaxis_cc*fac_r
        end if

        call nc_close(ncid)

    end subroutine vmecin

    !> Full-mesh lambda rows 1..kparb from the half-mesh rows with the map
    !> selected by lambda_half_to_full ('cubic' needs kparb >= 4, otherwise the
    !> two-point average with linear extrapolation at the edge). Row 0 (axis)
    !> is set to zero. The map is linear; boozer_vmec_ad tabulates its matrix.
    subroutine lambda_half_to_full_map(lhalf, lfull, kparb)
        use libneo_kinds, only: dp
        use new_vmec_stuff_mod, only: lambda_half_to_full
        integer, intent(in) :: kparb
        real(dp), intent(in) :: lhalf(:, 0:)
        real(dp), intent(inout) :: lfull(:, 0:)
        integer :: i

        lfull(:, 0) = 0d0
        if (trim(lambda_half_to_full) == 'cubic' .and. kparb >= 4) then
            call lambda_half_to_full_cubic(lhalf, lfull, kparb)
        else
            do i = 1, kparb - 1
                lfull(:, i) = 0.5d0*(lhalf(:, i + 1) + lhalf(:, i))
            end do
            lfull(:, kparb) = lhalf(:, kparb) + 0.5d0*(lhalf(:, kparb) - lhalf(:, kparb - 1))
        end if
    end subroutine lambda_half_to_full_map

    !> Full-mesh lambda at s_i = i*hs, i = 1..kparb, from the half-mesh rows
    !> 1..kparb at s = (j - 1/2)*hs by four-point Lagrange interpolation
    !> (centred inside, one-sided at the edge). Row 0 is left untouched.
    subroutine lambda_half_to_full_cubic(lhalf, lfull, kparb)
        use libneo_kinds, only: dp
        integer, intent(in) :: kparb
        real(dp), intent(in) :: lhalf(:, 0:)
        real(dp), intent(inout) :: lfull(:, 0:)

        integer :: i, j, k, jbeg
        real(dp) :: w(4)

        do i = 1, kparb
            jbeg = max(1, min(i - 1, kparb - 3))
            do j = 1, 4
                w(j) = 1d0
                do k = 1, 4
                    if (k /= j) w(j) = w(j)*(dble(i) - (dble(jbeg + k - 1) - 0.5d0)) &
                                       /(dble(jbeg + j - 1) - dble(jbeg + k - 1))
                end do
            end do
            lfull(:, i) = w(1)*lhalf(:, jbeg) + w(2)*lhalf(:, jbeg + 1) &
                          + w(3)*lhalf(:, jbeg + 2) + w(4)*lhalf(:, jbeg + 3)
        end do
    end subroutine lambda_half_to_full_cubic

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine stevvo(RT0, R0i, L1i, cbfi, bz0i, bf0)

        use libneo_kinds, only: dp
        use new_vmec_stuff_mod, only: nper, rmajor

        integer :: L1i
        real(dp) :: RT0, R0i, cbfi, bz0i, bf0

        associate (dummy => R0i)
        end associate
        associate (dummy => cbfi)
        end associate
        associate (dummy => bz0i)
        end associate
        associate (dummy => bf0)
        end associate

        L1i = nper
        RT0 = rmajor*1.d2

    end subroutine stevvo

end module vmecin_sub
