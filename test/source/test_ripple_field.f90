program test_ripple_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use analytical_tokamak_field
    implicit none

    type(analytical_circular_eq_t) :: eq_no_ripple, eq_with_ripple
    real(dp) :: R0, epsilon, A_param, B0
    real(dp) :: psi_target, dpsi_dR, dpsi_dZ, psi_val
    real(dp) :: R_flux, Z_scan, phi_scan
    real(dp) :: B_R, B_Z, B_phi_no_ripple, B_phi_ripple, B_mod
    integer :: i, j, nR, nZ, nphi
    integer :: unit_no_ripple, unit_ripple
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
    real(dp) :: R_axis, R_min, R_max, Z_min, Z_max
    real(dp) :: phi_min, phi_max
    logical :: found_flux_surface

    R0 = 6.2_dp
    epsilon = 0.32_dp
    A_param = -0.142_dp
    B0 = 5.3_dp

    call eq_no_ripple%init(R0, epsilon, A_in=A_param, B0_in=B0)

    call eq_with_ripple%init(R0, epsilon, A_in=A_param, B0_in=B0, &
                             Nripple_in=9, a0_in=2.0_dp, &
                             alpha0_in=2.0_dp, delta0_in=0.10_dp, z0_in=0.0_dp)

    R_axis = R0 + R0 * epsilon * 0.1_dp
    psi_target = eq_no_ripple%eval_psi(R_axis, 0.0_dp)

    nR = 100
    nZ = 100
    nphi = 72

    R_min = R0 - 2.5_dp
    R_max = R0 + 2.5_dp
    Z_min = -2.5_dp
    Z_max = 2.5_dp
    phi_min = 0.0_dp
    phi_max = 2.0_dp * pi

    found_flux_surface = .false.
    do i = 1, nR
        R_flux = R_min + (R_max - R_min) * real(i - 1, dp) / real(nR - 1, dp)
        psi_val = eq_no_ripple%eval_psi(R_flux, 0.0_dp)
        if (abs(psi_val - psi_target) < 1.0e-3_dp * abs(psi_target) .and. R_flux > R0) then
            found_flux_surface = .true.
            exit
        end if
    end do

    if (.not. found_flux_surface) then
        R_flux = R0 + 0.5_dp
        print *, "Warning: Could not find exact flux surface, using R =", R_flux
    else
        print *, "Found flux surface at R =", R_flux
    end if

    open(newunit=unit_no_ripple, file='bphi_no_ripple.csv', status='replace', action='write')
    write(unit_no_ripple, '(A)') 'phi,Z,Bphi'

    open(newunit=unit_ripple, file='bphi_with_ripple.csv', status='replace', action='write')
    write(unit_ripple, '(A)') 'phi,Z,Bphi'

    do j = 1, nZ
        Z_scan = Z_min + (Z_max - Z_min) * real(j - 1, dp) / real(nZ - 1, dp)
        do i = 1, nphi
            phi_scan = phi_min + (phi_max - phi_min) * real(i - 1, dp) / real(nphi - 1, dp)

            call eq_no_ripple%eval_bfield_ripple(R_flux, phi_scan, Z_scan, &
                                                 B_R, B_Z, B_phi_no_ripple, B_mod)
            write(unit_no_ripple, '(ES23.15E3,",",ES23.15E3,",",ES23.15E3)') &
                phi_scan, Z_scan, B_phi_no_ripple

            call eq_with_ripple%eval_bfield_ripple(R_flux, phi_scan, Z_scan, &
                                                   B_R, B_Z, B_phi_ripple, B_mod)
            write(unit_ripple, '(ES23.15E3,",",ES23.15E3,",",ES23.15E3)') &
                phi_scan, Z_scan, B_phi_ripple
        end do
    end do

    close(unit_no_ripple)
    close(unit_ripple)

    call eq_no_ripple%cleanup()
    call eq_with_ripple%cleanup()

    print *, "Data written to bphi_no_ripple.csv and bphi_with_ripple.csv"
    print *, "Ripple parameters: Nripple=9, delta0=0.10, alpha0=2.0"

end program test_ripple_field
