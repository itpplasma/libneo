!> Near-axis scalar-primary continuation through every straight-field-line route.
!> Tangents perturb m >= 2 only, so an O(1) profile tangent cannot hide a
!> continuation error. Cartesian contractions and shrinking-rho field variation
!> independently check the physical primal, alongside JVP/VJP and FD checks.
program test_boozer_axis_ad
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use boozer_hybrid, only: boozer_hybrid_t, hybrid_surface_t, &
                             boozer_hybrid_init, hybrid_surface, &
                             hybrid_surface_jvp, hybrid_surface_vjp, &
                             SFL_STRAIGHT_CURRENT, sfl_route_name
    implicit none
    real(dp), parameter :: TWOPI = 8.0_dp*atan(1.0_dp)
    real(dp), parameter :: STEPS(9) = [1.0e-2_dp, 1.0e-3_dp, 1.0e-4_dp, &
                                       1.0e-5_dp, 1.0e-6_dp, 3.0e-7_dp, &
                                       1.0e-7_dp, 1.0e-8_dp, 1.0e-9_dp]
    real(dp), parameter :: FD_WEIGHTS(4) = &
                           [4.0_dp/5, -1.0_dp/5, 4.0_dp/105, -1.0_dp/280]
    integer :: failures = 0

    call check_case('wout.nc', 'QA', .false.)
    call check_case('wout_li383.nc', 'li383', .true.)
    if (failures /= 0) error stop 'test_boozer_axis_ad FAILED'
    print '(a)', 'test_boozer_axis_ad PASSED'

contains

    subroutine check(label, error, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: error, tolerance
        print '(a, es11.3, a, es9.2)', trim(label)//': ', error, ' <= ', tolerance
        if (.not. ieee_is_finite(error) .or. error > tolerance) then
            failures = failures + 1
            print '(a)', 'FAIL '//trim(label)
        end if
    end subroutine check

    subroutine check_case(filename, label, pressure)
        character(*), intent(in) :: filename, label
        logical, intent(in) :: pressure
        type(boozer_hybrid_t) :: h, tangent, adjoint
        type(hybrid_surface_t) :: sv
        real(dp), allocatable :: yd(:), seed(:), fd(:), yp(:), ym(:)
        real(dp), allocatable :: sd(:, :), sb(:, :)
        real(dp) :: surfaces(3), fd_error(size(STEPS)), duality, source_error
        real(dp) :: components(4), component_best(4), best(4), best_step(4), err
        real(dp) :: io, ci, cg
        real(dp) :: sj, lhs, rhs, scale, eps
        integer :: route, j, k, i, l, nt, np, nout

        call boozer_hybrid_init(h, filename, 8, 6)
        nt = max(h%nt, 2*h%mnu + 8)
        np = max(h%np, 2*h%nnu + 8)
        if (h%axisym) np = 1
        surfaces = [0.25_dp, 0.75_dp, 1.25_dp]/real(h%ns - 1, dp)
        call check_regularity(h, label)
        nout = 6*h%mnboz + 3
        allocate (sd(6, h%mnboz), sb(6, h%mnboz), yd(nout), seed(nout), &
                  fd(nout), yp(nout), ym(nout))
        do route = 1, 5
            if (route == SFL_STRAIGHT_CURRENT .and. .not. pressure) cycle
            ! This source grid resolves both source and route Fourier modes.
            call boozer_hybrid_init(h, filename, 8, 6, n_theta=nt, n_phi=np, &
                                    route=route)
            tangent = h
            tangent%rmnc = 0.0_dp
            tangent%zmns = 0.0_dp
            tangent%lmns = 0.0_dp
            tangent%iotas = 0.0_dp
            do j = 0, h%ns - 1
                sj = real(j, dp)/real(h%ns - 1, dp)
                do i = 1, h%mnmax
                    if (abs(h%xm(i)) < 2) cycle
                    tangent%rmnc(i, j) = cos(real(i, dp))*(1 + 0.2_dp*sj)*h%rmnc(i, j)
                    tangent%zmns(i, j) = sin(real(i, dp))*(1 - 0.3_dp*sj)*h%zmns(i, j)
                    tangent%lmns(i, j) = cos(0.7_dp*i)*(1 + 0.4_dp*sj)*h%lmns(i, j)
                end do
            end do
            fd_error = 0.0_dp
            best = 0.0_dp
            duality = 0.0_dp
            source_error = 0.0_dp
            do j = 1, size(surfaces)
                call hybrid_surface(h, surfaces(j), sv)
                source_error = max(source_error, cartesian_error(h, sv))
                call hybrid_surface_jvp(h, sv, tangent%rmnc, tangent%zmns, &
                                        tangent%lmns, tangent%iotas, sd, io, ci, cg)
                yd = [reshape(sd, [6*h%mnboz]), io, ci, cg]
                if (sum(yd**2) <= tiny(1.0_dp)) &
                    error stop 'test_boozer_axis_ad: expected nonzero geometry tangent'
                component_best = huge(1.0_dp)
                best_step = 0.0_dp
                if (route == SFL_STRAIGHT_CURRENT) &
                    print '(a, i2, 4es14.5)', 'HamadaII surface, s/den/eta/|Jt|: ', &
                    j, surfaces(j), sv%den, maxval(abs(sv%eta)), sqrt(sum(yd**2))
                do i = 1, nout
                    seed(i) = sin(0.37_dp*i) + cos(0.61_dp*i)
                end do
                ! Give the three spectrum quantities comparable cotangents.
                sb = reshape(seed(:6*h%mnboz), shape(sb))
                do i = 1, 5, 2
                    scale = maxval(abs(sv%spec(i:i + 1, :)))
                    if (scale > 0.0_dp) sb(i:i + 1, :) = sb(i:i + 1, :)/scale
                end do
                seed(:6*h%mnboz) = reshape(sb, [6*h%mnboz])
                adjoint = h
                adjoint%rmnc = 0.0_dp
                adjoint%zmns = 0.0_dp
                adjoint%lmns = 0.0_dp
                adjoint%iotas = 0.0_dp
                call hybrid_surface_vjp(h, sv, sb, seed(nout - 2), seed(nout - 1), &
                                        seed(nout), adjoint%rmnc, adjoint%zmns, &
                                        adjoint%lmns, adjoint%iotas)
                lhs = dot_product(seed, yd)
                rhs = sum(adjoint%rmnc*tangent%rmnc) + sum(adjoint%zmns*tangent%zmns) &
                      + sum(adjoint%lmns*tangent%lmns)
                if (ieee_is_finite(rhs)) then
                    duality = max(duality, abs(lhs - rhs)/sqrt(sum(seed**2)*sum(yd**2)))
                else
                    duality = huge(1.0_dp)
                end if
                ! Eighth-order differences resolve the sensitive Hamada II map.
                ! Each physical output group and surface has its own roundoff
                ! optimum; every such minimum must meet the same strict bound.
                do k = 1, size(STEPS)
                    eps = STEPS(k)
                    fd = 0.0_dp
                    do l = 1, size(FD_WEIGHTS)
                        call perturbed(h, tangent, surfaces(j), l*eps, yp)
                        call perturbed(h, tangent, surfaces(j), -l*eps, ym)
                        fd = fd + FD_WEIGHTS(l)*(yp - ym)
                    end do
                    fd = fd/eps
                    err = derivative_error(fd, yd, h%mnboz, components)
                    where (components < component_best) best_step = eps
                    component_best = min(component_best, components)
                    fd_error(k) = max(fd_error(k), err)
                end do
                best = max(best, component_best)
                if (route == SFL_STRAIGHT_CURRENT) &
                    print '(a, 4es14.5)', 'HamadaII B/nu/grad/I,G best: ', &
                    component_best
                if (route == SFL_STRAIGHT_CURRENT) &
                    print '(a, 4es14.5)', 'HamadaII selected steps: ', best_step
            end do
            print '(a, 9es11.3)', &
                label//' '//trim(sfl_route_name(route))//' FD curve:', fd_error
            print '(a, 4es11.3)', &
                label//' '//trim(sfl_route_name(route))//' B/nu/grad/I,G minima:', best
            call check(label//' '//trim(sfl_route_name(route))//' JVP/FD', &
                       maxval(best), 1.0e-8_dp)
            call check(label//' '//trim(sfl_route_name(route))//' adjoint', &
                       duality, 1.0e-13_dp)
            call check(label//' '//trim(sfl_route_name(route))//' Cartesian B', &
                       source_error, 3.0e-13_dp)
        end do
    end subroutine check_case

    subroutine perturbed(h, tangent, s, step, y)
        type(boozer_hybrid_t), intent(in) :: h, tangent
        real(dp), intent(in) :: s, step
        real(dp), intent(out) :: y(:)
        type(boozer_hybrid_t) :: hp
        type(hybrid_surface_t) :: sv
        hp = h
        hp%rmnc = h%rmnc + step*tangent%rmnc
        hp%zmns = h%zmns + step*tangent%zmns
        hp%lmns = h%lmns + step*tangent%lmns
        call hybrid_surface(hp, s, sv)
        y = [reshape(sv%spec, [6*h%mnboz]), sv%iota, sv%ci, sv%cg]
    end subroutine perturbed

    real(dp) function derivative_error(fd, yd, nm, components) result(error)
        integer, intent(in) :: nm
        real(dp), intent(in) :: fd(:), yd(:)
        real(dp), intent(out) :: components(4)
        real(dp) :: delta(6, nm), tangent(6, nm), norm
        integer :: i
        if (.not. all(ieee_is_finite(fd)) .or. .not. all(ieee_is_finite(yd))) then
            error = huge(1.0_dp)
            components = error
            return
        end if
        delta = reshape(fd(:6*nm) - yd(:6*nm), [6, nm])
        tangent = reshape(yd(:6*nm), [6, nm])
        components = 0.0_dp
        do i = 1, 5, 2
            norm = sum(tangent(i:i + 1, :)**2)
            if (norm > tiny(1.0_dp)) then
                components((i + 1)/2) = sqrt(sum(delta(i:i + 1, :)**2)/norm)
            else if (sum(delta(i:i + 1, :)**2) > tiny(1.0_dp)) then
                components((i + 1)/2) = 1.0_dp
            end if
        end do
        norm = sum(yd(6*nm + 2:)**2)
        if (norm > tiny(1.0_dp)) then
            components(4) = sqrt(sum((fd(6*nm + 2:) - yd(6*nm + 2:))**2)/norm)
        else if (sum(fd(6*nm + 2:)**2) > tiny(1.0_dp)) then
            components(4) = 1.0_dp
        end if
        error = maxval(components)
    end function derivative_error

    real(dp) function cartesian_error(h, sv) result(error)
        type(boozer_hybrid_t), intent(in) :: h
        type(hybrid_surface_t), intent(in) :: sv
        real(dp) :: es(3), et(3), ep(3), normal(3), b(3), determinant, expected(4)
        real(dp) :: scale(4)
        integer :: p
        if (.not. all(ieee_is_finite(sv%f))) then
            error = huge(1.0_dp)
            return
        end if
        error = 0.0_dp
        scale = maxval(abs(sv%f), dim=2)
        do p = 1, size(sv%x, 2)
            ! Rotate the cylindrical basis to phi=0; dot products are invariant.
            es = [sv%x(2, p), 0.0_dp, sv%x(5, p)]
            et = [sv%x(3, p), 0.0_dp, sv%x(6, p)]
            ep = [sv%x(4, p), sv%x(1, p), sv%x(7, p)]
            normal = cross(et, ep)
            determinant = dot_product(es, normal)
            b = -h%phiedge/(TWOPI*determinant) &
                *((sv%iota - sv%x(9, p))*et + (1 + sv%x(8, p))*ep)
            expected = [sqrt(sum(b**2)), dot_product(b, et), dot_product(b, ep), &
                        sqrt(sum(normal**2))/abs(determinant)]
            error = max(error, maxval(abs(expected - sv%f(:, p))/scale))
        end do
    end function cartesian_error

    function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
             a(1)*b(2) - a(2)*b(1)]
    end function cross

    subroutine check_regularity(h, label)
        type(boozer_hybrid_t), intent(in) :: h
        character(*), intent(in) :: label
        type(hybrid_surface_t) :: sv
        real(dp), allocatable :: b(:, :)
        real(dp) :: variation(3), rho
        integer :: k, j
        allocate (b(h%nt, h%np))
        do k = 1, 3
            rho = 1.0e-3_dp/2.0_dp**(k - 1)
            call hybrid_surface(h, rho**2, sv, with_spectrum=.false.)
            b = reshape(sv%f(1, :), shape(b))
            do j = 1, h%np
                b(:, j) = b(:, j) - sum(b(:, j))/real(h%nt, dp)
            end do
            variation(k) = maxval(abs(b))/maxval(sv%f(1, :))
            call check(label//' near-axis Cartesian B', &
                       cartesian_error(h, sv), 3.0e-13_dp)
        end do
        print '(a, 3es11.3)', label//' poloidal |B| variation:', variation
        call check(label//' poloidal |B| variation contracts with rho', &
                   maxval(variation(2:)/variation(:2)), 0.60_dp)
    end subroutine check_regularity

end program test_boozer_axis_ad
