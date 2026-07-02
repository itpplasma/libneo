program test_ripple_field
    !> Quantitative checks for TF ripple perturbation in eval_bfield_ripple:
    !> the perturbation must be N-fold periodic in phi, have a peak-to-peak
    !> amplitude of the same order as delta0, and scale linearly with delta0.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use analytical_tokamak_field, only: analytical_circular_eq_t
    implicit none

    integer :: test_count, pass_count
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

    test_count = 0
    pass_count = 0

    print *, "Testing TF ripple perturbation (amplitude and periodicity)..."

    call test_ripple_periodicity()
    call test_ripple_amplitude_order_delta0()
    call test_ripple_amplitude_scales_with_delta0()
    call test_zero_ripple_recovers_axisymmetric()

    print *, ""
    print *, "===================="
    print *, "Test Summary:"
    print *, "  Total tests: ", test_count
    print *, "  Passed:      ", pass_count
    print *, "  Failed:      ", test_count - pass_count
    print *, "===================="

    if (pass_count /= test_count) then
        error stop "Some tests failed!"
    end if

contains

    subroutine init_default_ripple(eq, delta0)
        type(analytical_circular_eq_t), intent(out) :: eq
        real(dp), intent(in) :: delta0
        real(dp) :: R0, epsilon, A_param, B0

        R0 = 6.2_dp
        epsilon = 0.32_dp
        A_param = -0.142_dp
        B0 = 5.3_dp

        call eq%init(R0, epsilon, A_in=A_param, B0_in=B0, &
                     Nripple_in=9, a0_in=2.0_dp, alpha0_in=2.0_dp, &
                     delta0_in=delta0, z0_in=0.0_dp)
    end subroutine init_default_ripple

    subroutine test_ripple_periodicity()
        !> B_phi(R, phi, Z) must repeat with period 2*pi/Nripple, and must
        !> not be periodic with half that period (cos(N*phi) changes sign).
        type(analytical_circular_eq_t) :: eq
        real(dp) :: R_flux, Z_flux, period
        real(dp) :: B_R, B_Z, B_phi_a, B_phi_b, B_mod
        real(dp) :: phi
        integer :: i
        real(dp), parameter :: tol = 1.0e-10_dp

        call init_default_ripple(eq, 0.10_dp)

        R_flux = eq%R0 + 1.0_dp
        Z_flux = 0.3_dp
        period = 2.0_dp * pi / real(eq%Nripple, dp)

        do i = 0, 7
            phi = real(i, dp) * 0.37_dp

            call eq%eval_bfield_ripple(R_flux, phi, Z_flux, B_R, B_Z, B_phi_a, B_mod)
            call eq%eval_bfield_ripple(R_flux, phi + period, Z_flux, B_R, B_Z, B_phi_b, B_mod)
            call assert_close(B_phi_b, B_phi_a, tol, "B_phi periodic with 2*pi/Nripple")

            call eq%eval_bfield_ripple(R_flux, phi + 0.5_dp * period, Z_flux, B_R, B_Z, B_phi_b, B_mod)
            call assert_not_close(B_phi_b, B_phi_a, 1.0e-3_dp, &
                "B_phi NOT periodic with pi/Nripple (half period flips sign)")
        end do

        call eq%cleanup()
        print *, ""
    end subroutine test_ripple_periodicity

    subroutine test_ripple_amplitude_order_delta0()
        !> Peak-to-peak ripple modulation of B_phi, scanned over phi at a
        !> fixed (R, Z), should be the same order of magnitude as delta0
        !> relative to the axisymmetric B_phi (within a generous factor,
        !> since the divergence-free correction trades the naive
        !> "B_phi *= 1 + delta0*cos(N*phi)" model for an exact one).
        type(analytical_circular_eq_t) :: eq
        real(dp) :: R_flux, Z_flux, delta0
        real(dp) :: B_R, B_Z, B_phi, B_mod
        real(dp) :: phi, Bphi_min, Bphi_max, Bphi_axisym
        real(dp) :: peak_to_peak_fraction
        integer :: i, nphi

        call init_default_ripple(eq, 0.10_dp)
        delta0 = eq%delta0

        R_flux = eq%R0 + 1.0_dp
        Z_flux = 0.0_dp
        nphi = 360

        Bphi_min = huge(1.0_dp)
        Bphi_max = -huge(1.0_dp)
        do i = 0, nphi - 1
            phi = 2.0_dp * pi * real(i, dp) / real(nphi, dp)
            call eq%eval_bfield_ripple(R_flux, phi, Z_flux, B_R, B_Z, B_phi, B_mod)
            Bphi_min = min(Bphi_min, B_phi)
            Bphi_max = max(Bphi_max, B_phi)
        end do

        Bphi_axisym = eq%B0 * eq%R0 / R_flux
        peak_to_peak_fraction = (Bphi_max - Bphi_min) / Bphi_axisym

        test_count = test_count + 1
        if (peak_to_peak_fraction > 0.3_dp * delta0 .and. &
            peak_to_peak_fraction < 3.0_dp * delta0) then
            pass_count = pass_count + 1
        else
            print *, "  FAIL: peak-to-peak B_phi ripple fraction order delta0"
            print *, "    delta0:                  ", delta0
            print *, "    peak-to-peak fraction:    ", peak_to_peak_fraction
        end if

        call eq%cleanup()
        print *, ""
    end subroutine test_ripple_amplitude_order_delta0

    subroutine test_ripple_amplitude_scales_with_delta0()
        !> The ripple correction is linear in delta0: doubling delta0 must
        !> double the peak-to-peak B_phi modulation.
        type(analytical_circular_eq_t) :: eq_small, eq_large
        real(dp) :: R_flux, Z_flux
        real(dp) :: B_R, B_Z, B_phi, B_mod
        real(dp) :: phi, min_small, max_small, min_large, max_large
        real(dp) :: ratio
        integer :: i, nphi

        call init_default_ripple(eq_small, 0.05_dp)
        call init_default_ripple(eq_large, 0.10_dp)

        R_flux = eq_small%R0 + 1.0_dp
        Z_flux = 0.0_dp
        nphi = 360

        min_small = huge(1.0_dp); max_small = -huge(1.0_dp)
        min_large = huge(1.0_dp); max_large = -huge(1.0_dp)
        do i = 0, nphi - 1
            phi = 2.0_dp * pi * real(i, dp) / real(nphi, dp)

            call eq_small%eval_bfield_ripple(R_flux, phi, Z_flux, B_R, B_Z, B_phi, B_mod)
            min_small = min(min_small, B_phi)
            max_small = max(max_small, B_phi)

            call eq_large%eval_bfield_ripple(R_flux, phi, Z_flux, B_R, B_Z, B_phi, B_mod)
            min_large = min(min_large, B_phi)
            max_large = max(max_large, B_phi)
        end do

        ratio = (max_large - min_large) / (max_small - min_small)

        call assert_close(ratio, 2.0_dp, 1.0e-6_dp, &
            "doubling delta0 doubles peak-to-peak B_phi ripple")

        call eq_small%cleanup()
        call eq_large%cleanup()
        print *, ""
    end subroutine test_ripple_amplitude_scales_with_delta0

    subroutine test_zero_ripple_recovers_axisymmetric()
        !> Nripple=0 must reproduce the plain axisymmetric field exactly.
        type(analytical_circular_eq_t) :: eq
        real(dp) :: R0, epsilon, A_param, B0
        real(dp) :: R, phi, Z
        real(dp) :: B_R, B_Z, B_phi, B_mod
        real(dp) :: B_R_ax, B_Z_ax, B_phi_ax, B_mod_ax
        real(dp), parameter :: tol = 1.0e-12_dp

        R0 = 6.2_dp
        epsilon = 0.32_dp
        A_param = -0.142_dp
        B0 = 5.3_dp

        call eq%init(R0, epsilon, A_in=A_param, B0_in=B0)

        R = R0 + 1.0_dp
        phi = 1.23_dp
        Z = 0.4_dp

        call eq%eval_bfield_ripple(R, phi, Z, B_R, B_Z, B_phi, B_mod)
        call eq%eval_bfield(R, Z, B_R_ax, B_Z_ax, B_phi_ax, B_mod_ax)

        call assert_close(B_R, B_R_ax, tol, "Nripple=0: B_R matches axisymmetric")
        call assert_close(B_Z, B_Z_ax, tol, "Nripple=0: B_Z matches axisymmetric")
        call assert_close(B_phi, B_phi_ax, tol, "Nripple=0: B_phi matches axisymmetric")

        call eq%cleanup()
        print *, ""
    end subroutine test_zero_ripple_recovers_axisymmetric

    subroutine assert_close(actual, expected, tolerance, test_name)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: test_name
        real(dp) :: diff

        test_count = test_count + 1
        diff = abs(actual - expected)

        if (diff < tolerance) then
            pass_count = pass_count + 1
        else
            print *, "  FAIL: ", trim(test_name)
            print *, "    Expected: ", expected
            print *, "    Got:      ", actual
            print *, "    Diff:     ", diff
        end if
    end subroutine assert_close

    subroutine assert_not_close(actual, expected, tolerance, test_name)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: test_name
        real(dp) :: diff

        test_count = test_count + 1
        diff = abs(actual - expected)

        if (diff >= tolerance) then
            pass_count = pass_count + 1
        else
            print *, "  FAIL: ", trim(test_name)
            print *, "    Expected NOT close to: ", expected
            print *, "    Got:                   ", actual
            print *, "    Diff:                  ", diff
        end if
    end subroutine assert_not_close

end program test_ripple_field
