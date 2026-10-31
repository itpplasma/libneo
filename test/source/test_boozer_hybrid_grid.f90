!> Public R1 grid adapter: labels, units and radial covectors against the
!> Cartesian curl of the source vector potential and finite-difference maps.
program test_boozer_hybrid_grid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord, &
                          boozer_to_vmec, vmec_to_boozer, boozer_quality, &
                          BOOZER_LEGACY, sfl_quality_config_t, sfl_status_t
    use boozer_hybrid, only: boozer_hybrid_t, hybrid_surface_t, &
                             boozer_hybrid_init, hybrid_surface
    use boozer_coordinates_mod, only: use_B_r, use_del_tp_B, ns_B, hs_B, &
                                     n_theta_B, n_phi_B, h_theta_B, h_phi_B
    use new_vmec_stuff_mod, only: vmec_B_scale, vmec_RZ_scale
    implicit none

    real(dp), parameter :: PI = acos(-1.0_dp), TWOPI = 2.0_dp*PI
    real(dp), parameter :: SAMPLES(3) = [0.17_dp, 0.43_dp, 0.79_dp]
    real(dp), parameter :: DH = 1.0e-4_dp
    type(sfl_quality_config_t) :: coarse
    type(sfl_status_t) :: status
    real(dp) :: baseline(9, 3, 3), scaled(9, 3, 3), legacy(9)
    integer :: failures
    character(len=32) :: mode

    failures = 0
    use_B_r = .true.
    use_del_tp_B = .true.
    vmec_B_scale = 1.0_dp
    vmec_RZ_scale = 1.0_dp
    coarse%allow_poor_equilibrium = .true.
    coarse%override_reason = 'test_boozer_hybrid_grid: coarse public fixtures'
    if (command_argument_count() > 0) then
        call get_command_argument(1, mode)
        if (command_argument_count() /= 1 .or. mode /= '--offgrid-metric') &
            error stop 'usage: test_boozer_hybrid_grid.x [--offgrid-metric]'
        ! Larger spline tensors qualify off-grid metric accuracy on a cluster.
        call check_case('wout_li383.nc', baseline, refinement=10)
        call check_case('wout.nc', scaled, refinement=10)
        if (failures /= 0) error stop 'test_boozer_hybrid_grid FAILED'
        print '(a)', 'test_boozer_hybrid_grid: off-grid metric checks passed'
        stop
    end if
    call get_boozer_coordinates('wout_li383.nc', status=status)
    call expect('default quality gate rejects coarse source', status%code /= 0)
    call expect('failure is not overridden', .not. status%overridden)
    if (allocated(status%message)) then
        call expect('failure explains override', &
                    index(status%message, 'allow_poor_equilibrium') > 0)
    else
        call expect('failure carries a message', .false.)
    end if
    call get_boozer_coordinates('wout_li383.nc', construction=999, status=status)
    call expect('invalid construction is rejected', status%code /= 0)
    if (allocated(status%message)) then
        call expect('invalid construction has an actionable message', &
                    index(status%message, 'construction') > 0)
    else
        call expect('invalid construction carries a message', .false.)
    end if
    call check_case('wout_li383.nc', baseline)
    vmec_B_scale = 2.0_dp
    vmec_RZ_scale = 3.0_dp
    call check_case('wout_li383.nc', scaled)
    call check_scaling(baseline, scaled)
    vmec_B_scale = 1.0_dp
    vmec_RZ_scale = 1.0_dp
    call get_boozer_coordinates('wout_li383.nc', construction=BOOZER_LEGACY, &
                                status=status)
    call expect('explicit legacy remains available', status%code == 0)
    call sample(SAMPLES(2), 1.1_dp, 0.2_dp, legacy)
    call expect('explicit legacy returns finite fields', all(ieee_is_finite(legacy)))
    call check_case('wout_li383.nc', scaled, nodes=.true.)
    call check_case('wout.nc', baseline)
    call check_case('wout.nc', scaled, nodes=.true.)
    if (failures /= 0) error stop 'test_boozer_hybrid_grid FAILED'
    print '(a)', 'test_boozer_hybrid_grid: all checks passed'

contains

    subroutine expect(label, condition)
        character(*), intent(in) :: label
        logical, intent(in) :: condition
        if (.not. condition) then
            print '(a)', 'FAIL '//label
            failures = failures + 1
        end if
    end subroutine expect

    subroutine check(label, error, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: error, tolerance
        if (.not. (error <= tolerance)) then
            print '(a, 2es14.5)', 'FAIL '//label//': ', error, tolerance
            failures = failures + 1
        end if
    end subroutine check

    subroutine check_case(filename, values, nodes, refinement)
        character(*), intent(in) :: filename
        real(dp), intent(out) :: values(9, 3, 3)
        logical, intent(in), optional :: nodes
        integer, intent(in), optional :: refinement
        type(boozer_hybrid_t) :: h
        type(hybrid_surface_t) :: sv, shifted_surface(4)
        real(dp) :: s, rho, theta, phi, theta_b, phi_b, tb, pb, xyz(3), b(3)
        real(dp) :: shifted(3, 4), tangent(3), fine(3), dummy(3), grad_s
        real(dp) :: offsets(4), angle_error, field_error, radial_error, fd_error
        real(dp) :: label_error, flux, scale, iota, metric_error, signal
        integer :: is, ia, k, it, ip, p, refinement_used
        logical :: at_nodes

        refinement_used = 5
        if (present(refinement)) refinement_used = refinement
        at_nodes = .false.
        if (present(nodes)) at_nodes = nodes
        print '(a, i0, a, l1)', 'R1 grid adapter: '//filename//', refinement=', &
            refinement_used, ', grid nodes=', at_nodes
        values = 0.0_dp
        call get_boozer_coordinates(filename, quality=coarse, status=status, &
                                    grid_refinment=refinement_used)
        if (filename == 'wout_li383.nc') then
            call expect('explicit override retains failing code', status%code /= 0)
            call expect('explicit override permits construction', status%overridden)
            call expect('quality record retains override', boozer_quality%overridden)
            call expect('quality record retains reason', &
                        trim(boozer_quality%override_reason) == &
                        trim(coarse%override_reason))
        else
            call expect('quality gate permits construction', &
                        status%code == 0 .or. status%overridden)
        end if
        if (status%code /= 0 .and. .not. status%overridden) return
        call boozer_hybrid_init(h, filename, 8, 6)
        offsets = [-DH, DH, -DH/2, DH/2]
        angle_error = 0
        field_error = 0
        radial_error = 0
        fd_error = 0
        label_error = 0
        metric_error = 0
        signal = 0
        flux = -h%phiedge/TWOPI*1.0e8_dp*vmec_B_scale*vmec_RZ_scale**2
        do is = 1, size(SAMPLES)
            s = SAMPLES(is)
            if (at_nodes) then
                k = max(1, min(ns_B - 2, nint(sqrt(s)/hs_B)))
                s = (k*hs_B)**2
            end if
            rho = sqrt(s)
            call hybrid_surface(h, s, sv, with_spectrum=.false.)
            if (at_nodes) then
                do k = 1, 4
                    call hybrid_surface(h, (rho + offsets(k))**2, &
                                        shifted_surface(k), with_spectrum=.false.)
                end do
            end if
            do ia = 1, 3
                it = 1 + (ia*h%nt)/7
                ip = 1 + mod(ia*h%np/5, h%np)
                p = it + (ip - 1)*h%nt
                theta_b = sv%thb(p)
                phi_b = sv%phb(p)
                if (at_nodes) then
                    theta_b = ((ia*(n_theta_B - 1))/7)*h_theta_B
                    phi_b = ((ia*(n_phi_B - 1))/5)*h_phi_B
                end if
                call boozer_to_vmec(s, theta_b, phi_b, theta, phi)
                if (.not. at_nodes) then
                    angle_error = max(angle_error, &
                                      abs(wrap(theta - h%theta(it), TWOPI)), &
                                      abs(wrap(phi - h%phi(ip), TWOPI/h%nfp)))
                end if
                call vmec_to_boozer(s, theta, phi, tb, pb)
                call check('inverse map round trip', &
                           abs(wrap(tb - theta_b, TWOPI)) + &
                           abs(wrap(pb - phi_b, TWOPI/h%nfp)), 1.0e-10_dp)
                call source_at(h, sv%rv, sv%rd, theta, phi, xyz, b, grad_s)
                call sample(s, theta_b, phi_b, values(:, ia, is))
                scale = norm2(b)*1.0e4_dp*vmec_B_scale
                field_error = max(field_error, &
                                  abs(values(7, ia, is)/scale - 1.0_dp))
                iota = dot_product(sv%rv(:, 3), h%iotas)
                label_error = max(label_error, &
                                  abs(values(1, ia, is)/(flux*s) - 1), &
                                  abs(values(3, ia, is)/flux - 1), &
                                  abs(-values(4, ia, is)/values(3, ia, is) - iota))
                field_error = max(field_error, &
                    abs(values(5, ia, is) - sv%ci*1.0e6_dp &
                        *vmec_B_scale*vmec_RZ_scale)/(abs(sv%cg)*1.0e6_dp &
                        *vmec_B_scale*vmec_RZ_scale), &
                    abs(values(6, ia, is)/(sv%cg*1.0e6_dp &
                        *vmec_B_scale*vmec_RZ_scale) - 1))
                metric_error = max(metric_error, &
                    abs(values(9, ia, is)/(grad_s*0.01_dp/vmec_RZ_scale) - 1))
                if (at_nodes) then
                    ! The strict radial oracle uses the source inverse. The
                    ! compatibility rho spline has a separate derivative floor.
                    call source_inverse(h, sv, theta_b, phi_b, theta, phi)
                    call source_at(h, sv%rv, sv%rd, theta, phi, xyz, b, grad_s)
                    do k = 1, 4
                        call source_inverse(h, shifted_surface(k), theta_b, &
                                            phi_b, theta, phi)
                        call source_at(h, shifted_surface(k)%rv, &
                            shifted_surface(k)%rd, theta, phi, &
                            shifted(:, k), dummy, grad_s)
                    end do
                    tangent = (shifted(:, 2) - shifted(:, 1))/(2*DH)
                    fine = (shifted(:, 4) - shifted(:, 3))/DH
                    scale = norm2(b)*norm2(fine)*1.0e6_dp &
                            *vmec_B_scale*vmec_RZ_scale
                    fd_error = max(fd_error, norm2(fine - tangent)/norm2(fine))
                    radial_error = max(radial_error, &
                        abs(2*rho*values(8, ia, is) - dot_product(b, fine) &
                            *1.0e6_dp*vmec_B_scale*vmec_RZ_scale)/scale)
                    signal = max(signal, abs(dot_product(b, fine))*1.0e6_dp &
                                 *vmec_B_scale*vmec_RZ_scale/scale)
                end if
            end do
        end do
        print '(a, 4es12.4)', 'errors: angle, field, label, metric: ', &
            angle_error, field_error, label_error, metric_error
        call check('R1 source angle map', angle_error, 3.0e-5_dp)
        call check('Cartesian field and covariant flux functions', &
                   field_error, 3.0e-5_dp)
        call check('toroidal flux radial labels', label_error, 3.0e-6_dp)
        if (at_nodes .or. refinement_used == 10) &
            call check('source metric radial label', metric_error, 3.0e-5_dp)
        if (at_nodes) then
            print '(a, 2es12.4)', 'errors at grid nodes: source covector, FD: ', &
                radial_error, fd_error
            call check('source inverse-map finite-difference convergence', &
                       fd_error, 1.0e-6_dp)
            call check('Cartesian radial covector at grid nodes', &
                       radial_error, 3.0e-5_dp)
            call expect('radial oracle resolves a nonzero signal', signal > 1.0e-5_dp)
        end if
    end subroutine check_case

    !> Source geometry in SI. B follows curl(A), without a covariant-field or
    !> beta kernel. Rows are the source's radial interpolation contract.
    subroutine source_at(h, rv, rd, theta, phi, xyz, b, grad_s, angle_data)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: rv(0:, :), rd(0:, :), theta, phi
        real(dp), intent(out) :: xyz(3), b(3), grad_s
        real(dp), intent(out), optional :: angle_data(3)
        real(dp) :: r, z, rs, zs, rt, zt, rp, zp, lt, lp, c, sn, phase
        real(dp) :: rc, rr, zc, zz, lc, ll, rcs, rrs, zcs, zzs
        real(dp) :: es(3), et(3), ep(3), jacobian, psip, iota, lambda
        integer :: k, parity, full, half, m, n

        r = 0
        z = 0
        rs = 0
        zs = 0
        rt = 0
        zt = 0
        rp = 0
        zp = 0
        lt = 0
        lp = 0
        lambda = 0
        do k = 1, h%mnmax
            m = h%xm(k)
            n = h%xn(k)
            parity = mod(m, 2)
            full = 1 + parity
            half = 3 + parity
            phase = m*theta - n*phi
            c = cos(phase)
            sn = sin(phase)
            rc = dot_product(rv(:, full), h%rmnc(k, :))
            rr = dot_product(rv(:, full), h%rmns(k, :))
            zc = dot_product(rv(:, full), h%zmnc(k, :))
            zz = dot_product(rv(:, full), h%zmns(k, :))
            lc = dot_product(rv(:, half), h%lmnc(k, :))
            ll = dot_product(rv(:, half), h%lmns(k, :))
            rcs = dot_product(rd(:, full), h%rmnc(k, :))
            rrs = dot_product(rd(:, full), h%rmns(k, :))
            zcs = dot_product(rd(:, full), h%zmnc(k, :))
            zzs = dot_product(rd(:, full), h%zmns(k, :))
            r = r + rc*c + rr*sn
            z = z + zc*c + zz*sn
            rs = rs + rcs*c + rrs*sn
            zs = zs + zcs*c + zzs*sn
            rt = rt + m*(-rc*sn + rr*c)
            zt = zt + m*(-zc*sn + zz*c)
            rp = rp - n*(-rc*sn + rr*c)
            zp = zp - n*(-zc*sn + zz*c)
            lambda = lambda + lc*c + ll*sn
            lt = lt + m*(-lc*sn + ll*c)
            lp = lp - n*(-lc*sn + ll*c)
        end do
        xyz = [r*cos(phi), r*sin(phi), z]
        es = [rs*cos(phi), rs*sin(phi), zs]
        et = [rt*cos(phi), rt*sin(phi), zt]
        ep = [rp*cos(phi) - r*sin(phi), rp*sin(phi) + r*cos(phi), zp]
        jacobian = dot_product(es, cross(et, ep))
        psip = -h%phiedge/TWOPI
        iota = dot_product(rv(:, 3), h%iotas)
        b = psip*((iota - lp)*et + (1 + lt)*ep)/jacobian
        grad_s = norm2(cross(et, ep))/abs(jacobian)
        if (present(angle_data)) angle_data = [lambda, lt, lp]
    end subroutine source_at

    !> Invert the source R1 angles, independently of the compatibility spline.
    !> Cardinal interpolation resolves its periodic potential and derivatives.
    subroutine source_inverse(h, sv, theta_b, phi_b, theta, phi)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp), intent(in) :: theta_b, phi_b
        real(dp), intent(out) :: theta, phi
        real(dp) :: xyz(3), bvec(3), grad_s, map(3), nu, nut, nup
        real(dp) :: a, b, c, d, determinant, f1, f2
        integer :: iteration

        theta = theta_b
        phi = phi_b
        do iteration = 1, 30
            call source_at(h, sv%rv, sv%rd, theta, phi, xyz, bvec, grad_s, map)
            nu = source_value(h, sv%nu, theta, phi)
            nut = source_value(h, sv%nut, theta, phi)
            nup = source_value(h, sv%nup, theta, phi)
            f1 = theta + map(1) + sv%iota*nu - theta_b
            f2 = phi + nu - phi_b
            if (abs(f1) + abs(f2) < 1.0e-13_dp) return
            a = 1 + map(2) + sv%iota*nut
            b = map(3) + sv%iota*nup
            c = nut
            d = 1 + nup
            determinant = a*d - b*c
            theta = theta - (d*f1 - b*f2)/determinant
            phi = phi - (a*f2 - c*f1)/determinant
        end do
        call expect('exact source inverse converges', .false.)
    end subroutine source_inverse

    real(dp) function source_value(h, field, theta, phi) result(value)
        type(boozer_hybrid_t), intent(in) :: h
        real(dp), intent(in) :: field(:), theta, phi
        real(dp) :: wt(h%nt), wp(h%np), delta
        integer :: k

        wt = 1
        do k = 1, h%nt
            delta = wrap(theta - h%theta(k), TWOPI)
            if (abs(delta) > 1.0e-12_dp) wt(k) = &
                sin(h%nt*delta/2)/(h%nt*tan(delta/2))
        end do
        wp = 1
        if (h%np > 1) then
            do k = 1, h%np
                delta = wrap((phi - h%phi(k))*h%nfp, TWOPI)
                if (abs(delta) > 1.0e-12_dp) wp(k) = &
                    sin(h%np*delta/2)/(h%np*tan(delta/2))
            end do
        end if
        value = sum(matmul(wt, reshape(field, [h%nt, h%np]))*wp)
    end function source_value

    pure function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
             a(1)*b(2) - a(2)*b(1)]
    end function cross

    pure real(dp) function wrap(angle, period)
        real(dp), intent(in) :: angle, period
        wrap = modulo(angle + period/2, period) - period/2
    end function wrap

    subroutine sample(s, theta, phi, values)
        real(dp), intent(in) :: s, theta, phi
        real(dp), intent(out) :: values(9)
        real(dp) :: d2ap, d3ap, dit, d2it, dgt, d2gt, db(3), d2b(6)
        real(dp) :: dbr(3), d2br(6)
        call splint_boozer_coord(s, theta, phi, 0, values(1), values(2), &
             values(3), values(4), d2ap, d3ap, values(5), dit, d2it, &
             values(6), dgt, d2gt, values(7), db, d2b, values(8), dbr, d2br, &
             values(9))
    end subroutine sample

    subroutine check_scaling(base, changed)
        real(dp), intent(in) :: base(9, 3, 3), changed(9, 3, 3)
        real(dp), parameter :: factors(9) = [18.0_dp, 18.0_dp, 18.0_dp, &
                                             18.0_dp, 6.0_dp, 6.0_dp, &
                                             2.0_dp, 6.0_dp, 1.0_dp/3]
        real(dp) :: scale
        integer :: k
        do k = 1, 9
            scale = maxval(abs(base(k, :, :)))*factors(k)
            call check('VMEC field/length scaling', &
                       maxval(abs(changed(k, :, :) - factors(k)*base(k, :, :))) &
                       /max(scale, tiny(1.0_dp)), 1.0e-10_dp)
        end do
    end subroutine check_scaling
end program test_boozer_hybrid_grid
