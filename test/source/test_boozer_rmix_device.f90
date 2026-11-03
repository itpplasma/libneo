!> splint_boozer_rmix_device returns the radial-mixed field of
!> splint_boozer_coord_device (mode BOOZER_SECDERS_RADIAL_MIXED) through one
!> output record. Both run here on the host and must agree to round-off.
!> get_boozer_batch_splines exposes the splines that both evaluate.
program test_boozer_rmix_device
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
    use boozer_sub, only: get_boozer_coordinates, sync_boozer_state, &
                          splint_boozer_coord_device, splint_boozer_rmix_device, &
                          boozer_rmix_field_t, BOOZER_SECDERS_RADIAL_MIXED, &
                          BOOZER_LEGACY, get_boozer_batch_splines, boozer_state
    use batch_interpolate, only: BatchSplineData1D, BatchSplineData3D, &
                                 evaluate_batch_splines_1d_der2, &
                                 evaluate_batch_splines_3d
    implicit none

    real(dp), parameter :: tol = 1.0e-13_dp
    type(boozer_rmix_field_t) :: fld
    real(dp) :: Ath, Aph, dAth, dAph, d2Aph, d3Aph
    real(dp) :: Bt, dBt, d2Bt, Bp, dBp, d2Bp
    real(dp) :: Bm, dBm(3), d2Bm(6), Br, dBr(3), d2Br(6)
    real(dp) :: s, th, ph, err
    type(BatchSplineData1D), pointer :: aphi, bcov
    type(BatchSplineData3D), pointer :: field3d
    real(dp) :: y1(1), dy1(1), d2y1(1), y2(2), dy2(2), d2y2(2), y3(1)
    real(dp), parameter :: twopi = 2.0_dp*3.14159265358979_dp
    logical :: ready
    integer :: i

    call get_boozer_batch_splines(aphi, bcov, field3d, ready)
    if (ready .or. associated(aphi) .or. associated(bcov) .or. &
        associated(field3d)) then
        print *, 'FAIL: Boozer batch splines reported before construction'
        error stop 1
    end if

    use_B_r = .false.
    use_del_tp_B = .false.
    call get_boozer_coordinates("wout_li383.nc", construction=BOOZER_LEGACY)
    call sync_boozer_state
    call get_boozer_batch_splines(aphi, bcov, field3d, ready)
    if (.not. ready) then
        print *, 'FAIL: Boozer batch splines not reported after construction'
        error stop 1
    end if

    err = 0.0_dp
    do i = 1, 40
        s = 0.02_dp + 0.95_dp*real(i - 1, dp)/39.0_dp
        th = 0.37_dp*real(i, dp) - 3.0_dp
        ph = -0.53_dp*real(i, dp)
        call splint_boozer_coord_device(s, th, ph, BOOZER_SECDERS_RADIAL_MIXED, &
                                        Ath, Aph, dAth, dAph, d2Aph, d3Aph, &
                                        Bt, dBt, d2Bt, Bp, dBp, d2Bp, &
                                        Bm, dBm, d2Bm, Br, dBr, d2Br)
        call splint_boozer_rmix_device(s, th, ph, fld)
        call compare(fld%A_phi, Aph)
        call compare(fld%dA_phi_dr, dAph)
        call compare(fld%d2A_phi_dr2, d2Aph)
        call compare(fld%B_vartheta, Bt)
        call compare(fld%dB_vartheta, dBt)
        call compare(fld%d2B_vartheta, d2Bt)
        call compare(fld%B_varphi, Bp)
        call compare(fld%dB_varphi, dBp)
        call compare(fld%d2B_varphi, d2Bp)
        call compare(fld%Bmod, Bm)
        call compare(fld%dBmod(1), dBm(1))
        call compare(fld%dBmod(2), dBm(2))
        call compare(fld%dBmod(3), dBm(3))
        call compare(fld%d2Bmod(1), d2Bm(1))
        call compare(fld%d2Bmod(2), d2Bm(2))
        call compare(fld%d2Bmod(3), d2Bm(3))

        ! The exposed splines are the ones behind the field: A_phi in s,
        ! B_vartheta and B_varphi in rho = sqrt(s), |B| in (rho, angles).
        call evaluate_batch_splines_1d_der2(aphi, s, y1, dy1, d2y1)
        call compare(y1(1), Aph)
        call compare(dy1(1), dAph)
        call evaluate_batch_splines_1d_der2(bcov, sqrt(s), y2, dy2, d2y2)
        call compare(y2(1), Bt)
        call compare(y2(2), Bp)
        call evaluate_batch_splines_3d(field3d, [sqrt(s), modulo(th, twopi), &
                                       modulo(ph, twopi/real(boozer_state%nper, dp))], y3)
        call compare(y3(1), Bm)
    end do

    print '(A,ES10.3)', 'max relative difference = ', err
    if (err > tol) then
        print *, 'FAIL: splint_boozer_rmix_device differs from the radial-mixed mode'
        error stop 1
    end if
    print *, 'PASS'

contains

    subroutine compare(actual, expected)
        real(dp), intent(in) :: actual, expected

        err = max(err, abs(actual - expected)/max(abs(expected), tiny(1.0_dp)))
    end subroutine compare

end program test_boozer_rmix_device
