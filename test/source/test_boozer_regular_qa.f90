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
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord, &
                          BOOZER_HYBRID, BOOZER_LEGACY
    use boozer_hybrid, only: boozer_hybrid_t, hybrid_surface_t, &
                             boozer_hybrid_init, hybrid_surface
    use boozer_chartmap, only: export_boozer_chartmap, load_boozer_from_chartmap
    use boozer_regular, only: regular_field_t, set_boozer_regular_enabled, &
                              eval_boozer_regular, REGULAR_FIRST
    implicit none

    real(dp), parameter :: TWOPI = 6.283185307179586476925286766559_dp
    real(dp), parameter :: PHI0 = 0.3_dp
    !> Larmor-radius parameter of a fusion alpha in the converter's CGS units.
    real(dp), parameter :: RO0 = 3.0e5_dp
    integer, parameter :: N_THETA = 16
    integer, parameter :: CONSTRUCTIONS(2) = [BOOZER_HYBRID, BOOZER_LEGACY]
    character(6), parameter :: NAMES(2) = [character(6) :: 'hybrid', 'legacy']
    real(dp) :: physical_coefficient
    integer :: n_fail, construction, ic

    n_fail = 0
    use_B_r = .false.
    call set_boozer_regular_enabled(.true.)
    physical_coefficient = source_coefficient()
    print '(a,es18.10)', 'Cartesian source near-axis max|h_theta|/s: ', &
        physical_coefficient
    do ic = 1, size(CONSTRUCTIONS)
        construction = CONSTRUCTIONS(ic)
        call get_boozer_coordinates('wout.nc', radial_spline_order=5, &
                                    angular_spline_order=5, grid_refinment=5, &
                                    construction=construction)
        call check_source(trim(NAMES(ic))//' internal VMEC->Boozer')
        call export_boozer_chartmap('regular_qa.chartmap.nc')
        call load_boozer_from_chartmap('regular_qa.chartmap.nc')
        call check_source(trim(NAMES(ic))//' Boozer chartmap')
    end do
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
        real(dp) :: hth_ref, ratio_prev(2), ratio(2), pth, coefficient, previous
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
        previous = 0.0_dp
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
            coefficient = hmax(2)/s
            if (construction == BOOZER_HYBRID) then
                if (ir >= 4) call expect(label//': Cartesian h_theta coefficient', &
                    abs(coefficient/physical_coefficient - 1) <= 1.0e-2_dp)
                if (ir >= 7) call expect(label//': h_theta coefficient converges', &
                    abs(coefficient - previous) <= 1.0e-3_dp*abs(coefficient))
            else
                call expect(label//': regular h_theta O(s)', &
                            hmax(2) <= 2.0e2_dp*s*hth_ref + 1.0e-300_dp)
            end if
            previous = coefficient
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
        ! Both evaluations represent the same polar field. The hybrid regular
        ! representation resamples that field onto twice as many radial intervals.
        ! Against central differences of the (s, theta) evaluation on an 8x
        ! finer radial grid, d|B|/ds errs by up to 1.5e-2 (regular) and
        ! 1.7e-2 (legacy) of this scale; the radial interpolation error
        ! dominates, and the difference here falls to 6.6e-5 and 7.4e-6 on
        ! 2x and 4x finer radial grids.
        call expect(label//': grad |B| equals (s,theta) evaluation', &
                    err(6) <= 1.0e-2_dp)
    end subroutine check_against_legacy

    !> Independent Cartesian curl(A), with I=<B.e_theta> and the source angle
    !> map. Only the primary radial rows and periodic nu potential are shared.
    real(dp) function source_coefficient() result(coefficient)
        type(boozer_hybrid_t) :: h
        type(hybrid_surface_t) :: sv
        real(dp), allocatable :: amplitudes(:, :)
        real(dp) :: b(3), et(3), map(3), iota, ci, bmin, theta, phi
        real(dp) :: nu, nut, nup, a, bb, c, d, determinant, f1, f2
        real(dp), parameter :: S = 1.0e-8_dp
        integer :: k, m, full, half, it, ip, iteration

        call boozer_hybrid_init(h, 'wout.nc', 16, 16, n_theta=72, n_phi=72)
        if (h%lasym) error stop 'source_coefficient: symmetric QA fixture required'
        call hybrid_surface(h, S, sv, with_spectrum=.false.)
        allocate (amplitudes(5, h%mnmax))
        do k = 1, h%mnmax
            m = h%xm(k)
            full = 1 + mod(m, 2)
            half = 3 + mod(m, 2)
            if (m >= 2) then
                amplitudes(:, k) = [dot_product(sv%axis_rv(:, m, 1), h%rmnc(k, :)), &
                    dot_product(sv%axis_rv(:, m, 1), h%zmns(k, :)), &
                    dot_product(sv%axis_rv(:, m, 2), h%lmns(k, :)), &
                    dot_product(sv%axis_rd(:, m, 1), h%rmnc(k, :)), &
                    dot_product(sv%axis_rd(:, m, 1), h%zmns(k, :))]
            else
                amplitudes(:, k) = [dot_product(sv%rv(:, full), h%rmnc(k, :)), &
                    dot_product(sv%rv(:, full), h%zmns(k, :)), &
                    dot_product(sv%rv(:, half), h%lmns(k, :)), &
                    dot_product(sv%rd(:, full), h%rmnc(k, :)), &
                    dot_product(sv%rd(:, full), h%zmns(k, :))]
            end if
        end do
        iota = dot_product(sv%rv(:, 3), h%iotas)
        ci = 0
        do ip = 1, h%np
            do it = 1, h%nt
                call primary_field(h, amplitudes, iota, h%theta(it), h%phi(ip), &
                                   b, et, map)
                ci = ci + dot_product(b, et)
            end do
        end do
        ci = ci/real(h%nt*h%np, dp)
        bmin = huge(1.0_dp)
        do it = 0, N_THETA - 1
            theta = TWOPI*real(it, dp)/real(N_THETA, dp)
            phi = PHI0
            do iteration = 1, 30
                call primary_field(h, amplitudes, iota, theta, phi, b, et, map)
                nu = periodic_value(h, sv%nu, theta, phi)
                nut = periodic_value(h, sv%nut, theta, phi)
                nup = periodic_value(h, sv%nup, theta, phi)
                f1 = theta + map(1) + iota*nu - TWOPI*real(it, dp)/N_THETA
                f2 = phi + nu - PHI0
                if (abs(f1) + abs(f2) < 1.0e-13_dp) exit
                a = 1 + map(2) + iota*nut
                bb = map(3) + iota*nup
                c = nut
                d = 1 + nup
                determinant = a*d - bb*c
                theta = theta - (d*f1 - bb*f2)/determinant
                phi = phi - (a*f2 - c*f1)/determinant
            end do
            call expect('source Boozer inverse converges', iteration <= 30)
            call primary_field(h, amplitudes, iota, theta, phi, b, et, map)
            bmin = min(bmin, norm2(b))
        end do
        coefficient = 100*abs(ci)/(S*bmin) ! CGS covector/field converts m to cm.
        call expect('physical h_theta coefficient is positive', coefficient > 0)
    end function source_coefficient

    subroutine primary_field(h, amplitudes, iota, theta, phi, b, et, map)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: amplitudes(:, :), iota, theta, phi
        real(dp), intent(out) :: b(3), et(3), map(3)
        real(dp) :: ct(h%mnmax), sn(h%mnmax), es(3), ep(3)
        real(dp) :: r, rs, rt, rp, zs, zt, zp, lt, lp

        ct = cos(h%xm*theta - h%xn*phi)
        sn = sin(h%xm*theta - h%xn*phi)
        r = sum(amplitudes(1, :)*ct)
        rs = sum(amplitudes(4, :)*ct)
        rt = -sum(h%xm*amplitudes(1, :)*sn)
        rp = sum(h%xn*amplitudes(1, :)*sn)
        zs = sum(amplitudes(5, :)*sn)
        zt = sum(h%xm*amplitudes(2, :)*ct)
        zp = -sum(h%xn*amplitudes(2, :)*ct)
        lt = sum(h%xm*amplitudes(3, :)*ct)
        lp = -sum(h%xn*amplitudes(3, :)*ct)
        es = [rs*cos(phi), rs*sin(phi), zs]
        et = [rt*cos(phi), rt*sin(phi), zt]
        ep = [rp*cos(phi) - r*sin(phi), rp*sin(phi) + r*cos(phi), zp]
        b = -h%phiedge/(TWOPI*dot_product(es, cross(et, ep))) &
            *((iota - lp)*et + (1 + lt)*ep)
        map = [sum(amplitudes(3, :)*sn), lt, lp]
    end subroutine primary_field

    real(dp) function periodic_value(h, field, theta, phi) result(value)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: field(:), theta, phi
        real(dp) :: wt(h%nt), wp(h%np), delta
        integer :: k

        wt = 1
        wp = 1
        do k = 1, h%nt
            delta = modulo(theta - h%theta(k) + TWOPI/2, TWOPI) - TWOPI/2
            if (abs(delta) > 1.0e-12_dp) wt(k) = &
                sin(h%nt*delta/2)/(h%nt*tan(delta/2))
        end do
        do k = 1, h%np
            delta = modulo((phi - h%phi(k))*h%nfp + TWOPI/2, TWOPI) - TWOPI/2
            if (abs(delta) > 1.0e-12_dp) wp(k) = &
                sin(h%np*delta/2)/(h%np*tan(delta/2))
        end do
        value = sum(matmul(wt, reshape(field, [h%nt, h%np]))*wp)
    end function periodic_value

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
             a(1)*b(2) - a(2)*b(1)]
    end function cross

end program test_boozer_regular_qa
