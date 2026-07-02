module boozmn_reader
!> Load Boozer splines from a booz_xform boozmn NetCDF.
!>
!> Reads the harmonics via boozmn_file, evaluates Bmod on a uniform grid by
!> Fourier summation, and hands the result to build_boozer_from_chartmap so
!> the boozmn and chartmap paths share one spline backend.

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: load_boozer_from_boozmn

    real(dp), parameter :: TWOPI = 8.0_dp*atan(1.0_dp)
    real(dp), parameter :: GAUSS_CM2_PER_TM2 = 1.0e8_dp
    real(dp), parameter :: GAUSS_CM_PER_TM = 1.0e6_dp
    real(dp), parameter :: GAUSS_PER_T = 1.0e4_dp

contains

    subroutine load_boozer_from_boozmn(boozmn_filename, nrho_in, ntheta_in, nzeta_in)
        use boozmn_file, only: boozmn_data_t, read_boozmn
        use boozer_chartmap_types, only: boozer_chartmap_data_t
        use boozer_sub, only: build_boozer_from_chartmap

        character(len=*), intent(in) :: boozmn_filename
        integer, intent(in), optional :: nrho_in, ntheta_in, nzeta_in

        type(boozmn_data_t) :: b
        type(boozer_chartmap_data_t) :: d
        integer :: nrho, ntheta, nzeta
        integer :: ir, it, iz, mn, mn00
        integer :: k
        real(dp), allocatable :: rho_half(:), s_half(:)
        real(dp), allocatable :: theta(:), zeta(:), bmnc_out(:, :)
        real(dp), allocatable :: bmod_geom(:, :, :)
        real(dp) :: torflux_si, angle
        real(dp) :: rho_min, rho_max, s_min, s_max

        nrho = 30
        ntheta = 48
        nzeta = 96
        if (present(nrho_in)) nrho = nrho_in
        if (present(ntheta_in)) ntheta = ntheta_in
        if (present(nzeta_in)) nzeta = nzeta_in

        call read_boozmn(boozmn_filename, b)
        if (b%asym) then
            error stop "asymmetric (lasym=1) boozmn not yet supported"
        end if

        allocate (s_half(b%nsurf), rho_half(b%nsurf))
        do k = 1, b%nsurf
            s_half(k) = (real(b%jlist(k), dp) - 1.5_dp)/real(b%ns - 1, dp)
        end do
        rho_half = sqrt(s_half)

        torflux_si = -b%phi(b%ns)/TWOPI

        rho_min = 1.0e-3_dp
        rho_max = 1.0_dp
        s_min = rho_min**2
        s_max = rho_max**2

        d%n_rho = nrho
        d%n_s = nrho
        d%nfp = b%nfp
        d%torflux = torflux_si*GAUSS_CM2_PER_TM2
        d%rho_min = rho_min
        d%rho_max = rho_max
        d%h_rho = (rho_max - rho_min)/real(nrho - 1, dp)
        d%h_s = (s_max - s_min)/real(nrho - 1, dp)

        allocate (d%rho(nrho), d%s(nrho))
        allocate (d%A_phi(nrho), d%B_theta(nrho), d%B_phi(nrho))
        do ir = 1, nrho
            d%rho(ir) = rho_min + d%h_rho*real(ir - 1, dp)
            d%s(ir) = s_min + d%h_s*real(ir - 1, dp)
        end do

        call interp_1d(b%buco(b%jlist), rho_half, d%rho, b%nsurf, nrho, d%B_theta)
        call interp_1d(b%bvco(b%jlist), rho_half, d%rho, b%nsurf, nrho, d%B_phi)
        call iota_integral(b%iota(b%jlist), rho_half, d%s, b%nsurf, nrho, &
                           torflux_si, d%A_phi)
        d%B_theta = d%B_theta*GAUSS_CM_PER_TM
        d%B_phi = d%B_phi*GAUSS_CM_PER_TM
        d%A_phi = d%A_phi*GAUSS_CM2_PER_TM2

        allocate (theta(ntheta), zeta(nzeta), bmnc_out(nrho, b%nmodes))
        do it = 1, ntheta
            theta(it) = TWOPI*real(it - 1, dp)/real(ntheta, dp)
        end do
        do iz = 1, nzeta
            zeta(iz) = TWOPI/real(b%nfp, dp)*real(iz - 1, dp)/real(nzeta, dp)
        end do

        call interp_modes(b%bmnc, rho_half, d%rho, b%ixm, b%nmodes, b%nsurf, &
                          nrho, bmnc_out)

        allocate (bmod_geom(nrho, ntheta, nzeta))
        do ir = 1, nrho
            do it = 1, ntheta
                do iz = 1, nzeta
                    bmod_geom(ir, it, iz) = 0.0_dp
                    do mn = 1, b%nmodes
                        angle = real(b%ixm(mn), dp)*theta(it) &
                                - real(b%ixn(mn), dp)*zeta(iz)
                        bmod_geom(ir, it, iz) = bmod_geom(ir, it, iz) &
                                                + bmnc_out(ir, mn)*cos(angle)
                    end do
                    bmod_geom(ir, it, iz) = bmod_geom(ir, it, iz)*GAUSS_PER_T
                end do
            end do
        end do

        d%n_theta = ntheta + 1
        d%n_phi = nzeta + 1
        d%h_theta = theta(2) - theta(1)
        d%h_phi = zeta(2) - zeta(1)
        allocate (d%Bmod(nrho, d%n_theta, d%n_phi))
        d%Bmod(:, 1:ntheta, 1:nzeta) = bmod_geom
        d%Bmod(:, d%n_theta, 1:nzeta) = bmod_geom(:, 1, :)
        d%Bmod(:, 1:ntheta, d%n_phi) = bmod_geom(:, :, 1)
        d%Bmod(:, d%n_theta, d%n_phi) = bmod_geom(:, 1, 1)

        mn00 = 0
        do mn = 1, b%nmodes
            if (b%ixm(mn) == 0 .and. b%ixn(mn) == 0) then
                mn00 = mn
                exit
            end if
        end do
        if (mn00 > 0) then
            d%rmajor = b%rmnc(mn00, b%nsurf)
        else
            d%rmajor = 1.0_dp
        end if

        call build_boozer_from_chartmap(d)
    end subroutine load_boozer_from_boozmn

    !> Interpolate boozmn Fourier coefficients from the half grid to rho_out.
    !> bmnc_h is shaped (nmn, nsurf): Fortran order from NetCDF (comput_surfs, mn_mode).
    subroutine interp_modes(bmnc_h, rho_half, rho_out, ixm, nmn, nsurf, nrho_out, &
                            bmnc_out)
        integer, intent(in) :: nmn, nsurf, nrho_out
        real(dp), intent(in) :: bmnc_h(nmn, nsurf)
        real(dp), intent(in) :: rho_half(nsurf)
        real(dp), intent(in) :: rho_out(nrho_out)
        integer, intent(in) :: ixm(nmn)
        real(dp), intent(out) :: bmnc_out(nrho_out, nmn)

        integer :: ir, mn, k
        real(dp) :: rho, frac, ratio

        do ir = 1, nrho_out
            rho = rho_out(ir)
            if (rho <= rho_half(1)) then
                do mn = 1, nmn
                    if (rho_half(1) > 0.0_dp) then
                        ratio = rho/rho_half(1)
                        bmnc_out(ir, mn) = bmnc_h(mn, 1)*ratio**min(ixm(mn), 50)
                    else
                        bmnc_out(ir, mn) = bmnc_h(mn, 1)
                    end if
                end do
            else if (rho >= rho_half(nsurf)) then
                do mn = 1, nmn
                    bmnc_out(ir, mn) = bmnc_h(mn, nsurf)
                end do
            else
                k = 1
                do while (k < nsurf)
                    if (rho_half(k + 1) >= rho) exit
                    k = k + 1
                end do
                if (rho_half(k + 1) > rho_half(k)) then
                    frac = (rho - rho_half(k))/(rho_half(k + 1) - rho_half(k))
                else
                    frac = 0.0_dp
                end if
                do mn = 1, nmn
                    bmnc_out(ir, mn) = (1.0_dp - frac)*bmnc_h(mn, k) &
                                       + frac*bmnc_h(mn, k + 1)
                end do
            end if
        end do
    end subroutine interp_modes

    !> Interpolate a 1D radial profile from the half grid to rho_out.
    subroutine interp_1d(vals_h, rho_half, rho_out, nsurf, nrho_out, vals_out)
        integer, intent(in) :: nsurf, nrho_out
        real(dp), intent(in) :: vals_h(nsurf)
        real(dp), intent(in) :: rho_half(nsurf)
        real(dp), intent(in) :: rho_out(nrho_out)
        real(dp), intent(out) :: vals_out(nrho_out)

        integer :: ir, k
        real(dp) :: rho, frac

        do ir = 1, nrho_out
            rho = rho_out(ir)
            if (rho <= rho_half(1)) then
                vals_out(ir) = vals_h(1)
            else if (rho >= rho_half(nsurf)) then
                vals_out(ir) = vals_h(nsurf)
            else
                k = 1
                do while (k < nsurf)
                    if (rho_half(k + 1) >= rho) exit
                    k = k + 1
                end do
                if (rho_half(k + 1) > rho_half(k)) then
                    frac = (rho - rho_half(k))/(rho_half(k + 1) - rho_half(k))
                else
                    frac = 0.0_dp
                end if
                vals_out(ir) = (1.0_dp - frac)*vals_h(k) + frac*vals_h(k + 1)
            end if
        end do
    end subroutine interp_1d

    !> A_phi(s) = -torflux_si * integral_0^s iota(s') ds', in T*m^2.
    subroutine iota_integral(iota_h, rho_half, s_out, nsurf, nrho_out, torflux_si, &
                             A_phi_out)
        integer, intent(in) :: nsurf, nrho_out
        real(dp), intent(in) :: iota_h(nsurf)
        real(dp), intent(in) :: rho_half(nsurf)
        real(dp), intent(in) :: s_out(nrho_out)
        real(dp), intent(in) :: torflux_si
        real(dp), intent(out) :: A_phi_out(nrho_out)

        integer :: ir, k
        real(dp) :: iota_at_s, s, s_half_sq(nsurf), cum(nsurf)

        s_half_sq = rho_half**2

        cum(1) = 0.0_dp
        do k = 2, nsurf
            cum(k) = cum(k - 1) + 0.5_dp*(iota_h(k - 1) + iota_h(k)) &
                     *(s_half_sq(k) - s_half_sq(k - 1))
        end do

        do ir = 1, nrho_out
            s = s_out(ir)
            if (s <= s_half_sq(1)) then
                iota_at_s = iota_h(1)
                A_phi_out(ir) = -torflux_si*iota_at_s*s
            else if (s >= s_half_sq(nsurf)) then
                A_phi_out(ir) = -torflux_si*(cum(nsurf) + &
                                             iota_h(nsurf)*(s - s_half_sq(nsurf)))
            else
                k = 1
                do while (k < nsurf)
                    if (s_half_sq(k + 1) >= s) exit
                    k = k + 1
                end do
                if (s_half_sq(k + 1) > s_half_sq(k)) then
                    iota_at_s = iota_h(k) + (iota_h(k + 1) - iota_h(k))* &
                                (s - s_half_sq(k))/(s_half_sq(k + 1) - s_half_sq(k))
                else
                    iota_at_s = iota_h(k)
                end if
                A_phi_out(ir) = -torflux_si*(cum(k) + 0.5_dp*(iota_h(k) + iota_at_s)* &
                                             (s - s_half_sq(k)))
            end if
        end do
    end subroutine iota_integral

end module boozmn_reader
