!> Near-axis artifacts of the LandremanPaul2021 QA Boozer field, before
!> (splint_boozer_coord in (s, theta)) and after (eval_boozer_regular in
!> (X, Y)), for the internal VMEC->Boozer transform and for a Boozer chartmap
!> exported from it.
!>
!> Printed per s: the poloidal spread of |B| over 16 angles, max |h_theta|,
!> the m = 0 cone coefficient (theta-averaged d|B|/drho)/B0, and p_theta/s
!> for v_par = +-1. Asserted for the regular field: |B| spread O(rho),
!> h_theta and A_theta O(s), vanishing cone, p_theta/s finite with one sign,
!> and agreement with the (s, theta) evaluation away from the axis.
program test_boozer_regular_qa
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_coordinates_mod, only: use_B_r
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord
    use boozer_chartmap, only: export_boozer_chartmap, load_boozer_from_chartmap
    use boozer_regular, only: regular_field_t, set_boozer_regular_enabled, &
                              eval_boozer_regular, REGULAR_FIRST
    implicit none

    real(dp), parameter :: TWOPI = 6.283185307179586476925286766559_dp
    real(dp), parameter :: PHI0 = 0.3_dp
    !> Larmor-radius parameter of a fusion alpha in the converter's CGS units.
    real(dp), parameter :: RO0 = 3.0e5_dp
    integer, parameter :: N_THETA = 16
    integer :: n_fail

    n_fail = 0
    use_B_r = .false.
    call set_boozer_regular_enabled(.true.)
    call get_boozer_coordinates('wout.nc', radial_spline_order=5, &
                                angular_spline_order=5, grid_refinment=5)
    call check_source('internal VMEC->Boozer')
    call export_boozer_chartmap('regular_qa.chartmap.nc')
    call load_boozer_from_chartmap('regular_qa.chartmap.nc')
    call check_source('Boozer chartmap')
    if (n_fail > 0) then
        print *, 'test_boozer_regular_qa: FAILED checks:', n_fail
        error stop 1
    end if
    print *, 'test_boozer_regular_qa: all checks passed'

contains

    subroutine expect(label, ok)
        character(*), intent(in) :: label
        logical, intent(in) :: ok

        if (.not. ok) then
            print '(2a)', 'FAIL ', label
            n_fail = n_fail + 1
        end if
    end subroutine expect

    !> Legacy evaluation in (s, theta, phi): |B|, d|B|/d(s, theta, phi),
    !> h_theta, h_phi, A_theta, A_phi.
    subroutine legacy(s, theta, phi, b, db, hth, hph, ath, aph)
        real(dp), intent(in) :: s, theta, phi
        real(dp), intent(out) :: b, db(3), hth, hph, ath, aph
        real(dp) :: dath, daph, d2aph, d3aph, bth, dbth, d2bth, bph, dbph, d2bph
        real(dp) :: d2b(6), br, dbr(3), d2br(6)

        call splint_boozer_coord(s, theta, phi, 0, ath, aph, dath, daph, d2aph, &
                                 d3aph, bth, dbth, d2bth, bph, dbph, d2bph, b, db, &
                                 d2b, br, dbr, d2br)
        hth = bth/b
        hph = bph/b
    end subroutine legacy

    subroutine check_source(label)
        character(*), intent(in) :: label
        type(regular_field_t) :: f, f_axis
        real(dp) :: s, rho, theta, b0, grad0, b, db(3), hth, hph, ath, aph
        real(dp) :: spread(2), hmax(2), cone(2), bmin(2), bmax(2), pmin(2), pmax(2)
        real(dp) :: hth_ref, ratio_prev(2), ratio(2), pth
        integer :: ir, it, iv

        print '(/,2a)', 'source: ', label
        call eval_boozer_regular(0.0_dp, 0.0_dp, PHI0, REGULAR_FIRST, f_axis)
        b0 = f_axis%Bmod
        grad0 = hypot(f_axis%dBmod(1), f_axis%dBmod(2))
        call eval_boozer_regular(0.1_dp, 0.0_dp, PHI0, REGULAR_FIRST, f)
        hth_ref = abs(f%hth)
        print '(a)', '       s  spread(B)/B0     max|hth|  cone/B0      '// &
            'p_th/s range (v_par=+-1)  [legacy | regular]'
        ratio_prev = 0.0_dp
        do ir = 2, 10
            s = 10.0_dp**(-ir)
            rho = sqrt(s)
            bmin = huge(1.0_dp)
            bmax = -huge(1.0_dp)
            hmax = 0.0_dp
            cone = 0.0_dp
            pmin = huge(1.0_dp)
            pmax = -huge(1.0_dp)
            do it = 0, N_THETA - 1
                theta = TWOPI*real(it, dp)/real(N_THETA, dp)
                call legacy(s, theta, PHI0, b, db, hth, hph, ath, aph)
                bmin(1) = min(bmin(1), b)
                bmax(1) = max(bmax(1), b)
                hmax(1) = max(hmax(1), abs(hth))
                cone(1) = cone(1) + 2.0_dp*rho*db(1)/real(N_THETA, dp)
                do iv = -1, 1, 2
                    pth = (hth*real(iv, dp) + ath/RO0)/s
                    pmin(1) = min(pmin(1), pth)
                    pmax(1) = max(pmax(1), pth)
                end do
                call eval_boozer_regular(rho*cos(theta), rho*sin(theta), PHI0, &
                                         REGULAR_FIRST, f)
                bmin(2) = min(bmin(2), f%Bmod)
                bmax(2) = max(bmax(2), f%Bmod)
                hmax(2) = max(hmax(2), abs(f%hth))
                cone(2) = cone(2) + (cos(theta)*f%dBmod(1) &
                                     + sin(theta)*f%dBmod(2))/real(N_THETA, dp)
                do iv = -1, 1, 2
                    pth = (f%hth*real(iv, dp) + f%Ath/RO0)/s
                    pmin(2) = min(pmin(2), pth)
                    pmax(2) = max(pmax(2), pth)
                end do
            end do
            spread = (bmax - bmin)/b0
            cone = cone/b0
            print '(es8.1,2(es10.2,es10.2,es10.2,2es11.3,a))', s, &
                spread(1), hmax(1), cone(1), pmin(1), pmax(1), ' |', &
                spread(2), hmax(2), cone(2), pmin(2), pmax(2), ''
            call expect(label//': regular |B| spread O(rho)', &
                        spread(2) <= 2.5_dp*rho*grad0/b0 + 1.0e-14_dp)
            call expect(label//': regular h_theta O(s)', &
                        hmax(2) <= 2.0e2_dp*s*hth_ref + 1.0e-300_dp)
            call expect(label//': regular cone vanishes', abs(cone(2)) <= 1.0e1_dp*rho)
            call expect(label//': regular p_theta/s one sign', &
                        pmin(2)*pmax(2) > 0.0_dp)
            ratio = [pmin(2), pmax(2)]
            if (ir > 6) call expect(label//': regular p_theta/s converges', &
                                    all(abs(ratio - ratio_prev) <= &
                                        1.0e-3_dp*abs(ratio)))
            ratio_prev = ratio
        end do
        call expect(label//': regular axis A_theta and h_theta vanish', &
                    f_axis%Ath == 0.0_dp .and. f_axis%hth == 0.0_dp)
        call check_against_legacy(label)
    end subroutine check_source

    !> Agreement with the (s, theta) evaluation away from the axis; derivatives
    !> are compared through the chain rule d/ds = (X d/dX + Y d/dY)/(2 s) and
    !> d/dtheta = X d/dY - Y d/dX.
    subroutine check_against_legacy(label)
        character(*), intent(in) :: label
        type(regular_field_t) :: f
        real(dp) :: s, theta, phi, x, y, b, db(3), hth, hph, ath, aph
        real(dp) :: err(6), dreg(3), scale_db(3)
        integer :: k

        err = 0.0_dp
        scale_db = 0.0_dp
        do k = 1, 64
            s = 0.04_dp + 0.85_dp*real(mod(7*k, 64), dp)/64.0_dp
            theta = 0.37_dp*real(k, dp)
            phi = 0.53_dp*real(k, dp)
            x = sqrt(s)*cos(theta)
            y = sqrt(s)*sin(theta)
            call legacy(s, theta, phi, b, db, hth, hph, ath, aph)
            call eval_boozer_regular(x, y, phi, REGULAR_FIRST, f)
            dreg(1) = (x*f%dBmod(1) + y*f%dBmod(2))/(2.0_dp*s)
            dreg(2) = x*f%dBmod(2) - y*f%dBmod(1)
            dreg(3) = f%dBmod(3)
            err(1) = max(err(1), abs(f%Bmod/b - 1.0_dp))
            err(2) = max(err(2), abs(f%hth - hth)/abs(hph))
            err(3) = max(err(3), abs(f%hph/hph - 1.0_dp))
            err(4) = max(err(4), abs(f%Ath - ath)/abs(ath))
            err(5) = max(err(5), abs(f%Aph - aph)/max(abs(aph), abs(ath)))
            err(6) = max(err(6), maxval(abs(dreg - db)))
            scale_db = max(scale_db, abs(db))
        end do
        err(6) = err(6)/maxval(scale_db)
        print '(a,6es10.2)', 'max rel. difference to (s,theta) evaluation, '// &
            's in [0.04,0.89]: B, hth, hph, Ath, Aph, dB:', err
        call expect(label//': |B| equals (s,theta) evaluation', err(1) <= 1.0e-5_dp)
        call expect(label//': h equals (s,theta) evaluation', &
                    err(2) <= 1.0e-5_dp .and. err(3) <= 1.0e-5_dp)
        ! A_phi(s) is resampled onto the rho nodes as A0 + s a(rho).
        call expect(label//': A equals (s,theta) evaluation', &
                    err(4) <= 1.0e-12_dp .and. err(5) <= 1.0e-8_dp)
        ! Both evaluations interpolate the same nodal |B| on 50 radial nodes.
        ! Against central differences of the (s, theta) evaluation on an 8x
        ! finer radial grid, d|B|/ds errs by up to 1.5e-2 (regular) and
        ! 1.7e-2 (legacy) of this scale; the radial interpolation error
        ! dominates, and the difference here falls to 6.6e-5 and 7.4e-6 on
        ! 2x and 4x finer radial grids.
        call expect(label//': grad |B| equals (s,theta) evaluation', &
                    err(6) <= 1.0e-2_dp)
    end subroutine check_against_legacy

end program test_boozer_regular_qa
