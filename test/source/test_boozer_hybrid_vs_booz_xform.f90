!> Compare the hybrid Boozer spline field of get_boozer_coordinates with an
!> independent booz_xform transform of the same wout on outer half-grid
!> surfaces of the Landreman-Paul QH equilibrium (nfp = 4).
!>
!> The consumer spline grid is filled from a hybrid transform whose output
!> cutoff is about 4 mpol. Its quadrature must resolve that cutoff, or the top
!> modes alias. They are small in |B| but enter dB/dphi with the factor n nfp.
!> A quadrature of 2 mboz + 8 points gives 1.4e-2 to 2.2e-2 in dB/dphi here.
program test_boozer_hybrid_vs_booz_xform
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_coordinates_mod, only: use_B_r
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord, &
                          sfl_quality_config_t
    use boozmn_file, only: boozmn_data_t, read_boozmn
    implicit none

    character(len=*), parameter :: wout_file = 'wout_qh.nc'
    character(len=*), parameter :: boozmn_file_name = 'boozmn_qh.nc'
    real(dp), parameter :: TWOPI = 8.0_dp*atan(1.0_dp)
    real(dp), parameter :: GAUSS_PER_T = 1.0e4_dp
    !> Relative to the surface maximum of each quantity. With the resolved
    !> quadrature the differences are at most 6e-5 in |B| and 3.1e-3 in the
    !> angular derivatives, unchanged from grid_refinement 3 to 6 and from
    !> mboz = 32 to 48.
    real(dp), parameter :: TOL_B = 2.0e-4_dp, TOL_DB = 5.0e-3_dp
    integer, parameter :: NTH = 8, NPH = 8

    type(boozmn_data_t) :: d
    type(sfl_quality_config_t) :: quality
    integer :: k, failures

    call read_boozmn(boozmn_file_name, d)
    ! The ns = 50 fixture misses the default eps_fb threshold (1.04e-3 against
    ! 1.0e-3). Both transforms here read the same wout, so the comparison does
    ! not depend on how well that equilibrium is converged.
    quality%allow_poor_equilibrium = .true.
    quality%override_reason = 'test_boozer_hybrid_vs_booz_xform: same-wout '// &
                              'comparison with booz_xform'
    use_B_r = .true.
    call get_boozer_coordinates(wout_file, 5, 5, 3, quality=quality)

    failures = 0
    do k = 1, d%nsurf
        call check_surface(k)
    end do
    if (failures /= 0) error stop 'test_boozer_hybrid_vs_booz_xform FAILED'
    print '(a)', 'test_boozer_hybrid_vs_booz_xform: all checks passed'

contains

    subroutine check_surface(k)
        integer, intent(in) :: k
        real(dp) :: s, theta, phi, arg, ref(3), hyb(3), err(3), scale(3)
        real(dp) :: A_theta, A_phi, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, &
                    d3A_phi_dr3, B_vt, dB_vt, d2B_vt, B_vp, dB_vp, d2B_vp, &
                    Bmod, B_r, sqrt_g_ss
        real(dp) :: dBmod(3), dB_r(3), d2Bmod(6), d2B_r(6)
        integer :: i, j, mn

        s = (real(d%jlist(k), dp) - 1.5_dp)/real(d%ns - 1, dp)
        err = 0.0_dp
        scale = 0.0_dp
        do j = 0, NPH - 1
            phi = TWOPI*(real(j, dp) + 0.37_dp)/real(NPH*d%nfp, dp)
            do i = 0, NTH - 1
                theta = TWOPI*(real(i, dp) + 0.21_dp)/real(NTH, dp)
                ref = 0.0_dp
                do mn = 1, d%nmodes
                    arg = d%ixm(mn)*theta - d%ixn(mn)*phi
                    ref(1) = ref(1) + d%bmnc(mn, k)*cos(arg)
                    ref(2) = ref(2) - d%ixm(mn)*d%bmnc(mn, k)*sin(arg)
                    ref(3) = ref(3) + d%ixn(mn)*d%bmnc(mn, k)*sin(arg)
                end do
                call splint_boozer_coord(s, theta, phi, 0, A_theta, A_phi, &
                                         dA_theta_dr, dA_phi_dr, d2A_phi_dr2, &
                                         d3A_phi_dr3, B_vt, dB_vt, d2B_vt, &
                                         B_vp, dB_vp, d2B_vp, Bmod, dBmod, &
                                         d2Bmod, B_r, dB_r, d2B_r, sqrt_g_ss)
                hyb = [Bmod, dBmod(2), dBmod(3)]/GAUSS_PER_T
                err = max(err, abs(hyb - ref))
                scale = max(scale, abs(ref))
            end do
        end do
        err = err/scale
        print '(a,f7.4,a,3es11.3)', 's = ', s, &
            '  max rel. diff. |B|, dB/dtheta, dB/dphi:', err
        if (err(1) > TOL_B .or. any(err(2:3) > TOL_DB)) then
            print '(a,f7.4)', 'FAIL: hybrid field differs from booz_xform at s = ', s
            failures = failures + 1
        end if
    end subroutine check_surface

end program test_boozer_hybrid_vs_booz_xform
