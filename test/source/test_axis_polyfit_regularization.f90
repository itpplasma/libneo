program test_axis_polyfit_regularization
    !> rho**|m| * polyfit(s) near-axis continuation (s_to_rho_polyfit): a harmonic
    !> with m>0 vanishes at the axis, an exact Zernike radial input c = rho**m *
    !> P(s) is reproduced (the reduced amplitude is a low-degree polynomial in s),
    !> the continuation is continuous at the anchor, and unlike the pure power-law
    !> path it does not amplify noise into near-axis oscillations.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use new_vmec_stuff_mod, only: ns_s, s_axis_heal, rho_axis_heal
    use spline_vmec_sub, only: s_to_rho_polyfit, axis_anchor_index
    implicit none

    integer, parameter :: ns = 100, nrho = 100, ndeg = 3
    real(dp), parameter :: amp = 0.75_dp

    ns_s = 5
    s_axis_heal = 1.0e-2_dp
    rho_axis_heal = -1.0_dp

    call test_axis_vanishing()
    call test_zernike_reproduction()
    call test_anchor_continuity()
    call test_no_noise_oscillation()
    call test_axis_anchor_pinning()

    print *, 'axis polyfit regularization: all checks passed'

contains

    subroutine test_axis_vanishing()
        integer :: i_anchor
        real(dp) :: arr_in(ns), arr_out(nrho)
        i_anchor = axis_anchor_index(ns)
        call fill_zernike(2, arr_in)
        call s_to_rho_polyfit(2, ns, nrho, i_anchor, ndeg, arr_in, arr_out)
        if (abs(arr_out(1)) > 1.0e-13_dp) then
            print *, 'arr_out at rho=0 =', arr_out(1)
            error stop 'm>0 harmonic does not vanish at the axis'
        end if
    end subroutine test_axis_vanishing

    !> c(rho) = amp*rho**m*(1 + 0.5 s + 0.2 s**2) has reduced amplitude
    !> c/rho**m = amp*(1 + 0.5 s + 0.2 s**2), a degree-2 polynomial in s, so a
    !> degree-3 polyfit must reproduce it on the whole rho grid.
    subroutine test_zernike_reproduction()
        integer, parameter :: ms(3) = [1, 2, 4]
        integer :: im, m, i_anchor, irho
        real(dp) :: arr_in(ns), arr_out(nrho)
        real(dp) :: rho, s, expected, max_err
        i_anchor = axis_anchor_index(ns)
        do im = 1, size(ms)
            m = ms(im)
            call fill_zernike(m, arr_in)
            call s_to_rho_polyfit(m, ns, nrho, i_anchor, ndeg, arr_in, arr_out)
            max_err = 0.0_dp
            do irho = 1, nrho
                rho = real(irho - 1, dp)/real(nrho - 1, dp)
                s = rho*rho
                expected = amp*rho**m*(1.0_dp + 0.5_dp*s + 0.2_dp*s*s)
                max_err = max(max_err, abs(arr_out(irho) - expected))
            end do
            if (max_err > 1.0e-8_dp) then
                print *, 'm =', m, ' max |arr_out - zernike| =', max_err
                error stop 'Zernike radial polynomial not reproduced'
            end if
        end do
    end subroutine test_zernike_reproduction

    subroutine test_anchor_continuity()
        integer :: i_anchor, irho, irho_below, irho_above
        real(dp) :: arr_in(ns), arr_out(nrho), rho, rho_anchor, jump
        i_anchor = axis_anchor_index(ns)
        rho_anchor = sqrt(real(i_anchor - 1, dp)/real(ns - 1, dp))
        call fill_zernike(2, arr_in)
        call s_to_rho_polyfit(2, ns, nrho, i_anchor, ndeg, arr_in, arr_out)
        irho_below = 0; irho_above = 0
        do irho = 1, nrho
            rho = real(irho - 1, dp)/real(nrho - 1, dp)
            if (rho < rho_anchor) irho_below = irho
            if (rho >= rho_anchor .and. irho_above == 0) irho_above = irho
        end do
        jump = abs(arr_out(irho_above) - arr_out(irho_below))
        if (jump > 5.0e-3_dp) then
            print *, 'jump across anchor =', jump
            error stop 'discontinuity in polyfit continuation at the anchor'
        end if
    end subroutine test_anchor_continuity

    !> Reduced-amplitude noise on the reliable surfaces must not blow up into
    !> a near-axis oscillation: the least-squares fit smooths it. The second
    !> difference of arr_out below the anchor stays comparable to the input
    !> perturbation, not orders of magnitude larger (the power-law failure mode).
    subroutine test_no_noise_oscillation()
        integer :: i_anchor, irho, n_below
        real(dp) :: arr_in(ns), arr_out(nrho)
        real(dp) :: rho, rho_anchor, d2, max_d2
        i_anchor = axis_anchor_index(ns)
        rho_anchor = sqrt(real(i_anchor - 1, dp)/real(ns - 1, dp))
        call fill_zernike(1, arr_in)
        ! add a small deterministic perturbation to the reliable surfaces
        do irho = i_anchor, ns
            arr_in(irho) = arr_in(irho) + 1.0e-3_dp*amp*sin(7.0_dp*real(irho, dp))
        end do
        call s_to_rho_polyfit(1, ns, nrho, i_anchor, ndeg, arr_in, arr_out)
        max_d2 = 0.0_dp
        n_below = 0
        do irho = 2, nrho - 1
            rho = real(irho - 1, dp)/real(nrho - 1, dp)
            if (rho >= rho_anchor) cycle
            d2 = abs(arr_out(irho + 1) - 2.0_dp*arr_out(irho) + arr_out(irho - 1))
            max_d2 = max(max_d2, d2)
            n_below = n_below + 1
        end do
        if (n_below > 0 .and. max_d2 > 1.0e-2_dp) then
            print *, 'max |second difference| below anchor =', max_d2
            error stop 'polyfit continuation oscillates near the axis'
        end if
    end subroutine test_no_noise_oscillation

    !> The m=0 axis pin: with an explicit anchor, s_to_rho_polyfit forces the
    !> healed amplitude at rho=0 to equal the supplied axis value (raxis_cc /
    !> zaxis_cs in production), regardless of what the bare polyfit extrapolates.
    !> For m/=0 the anchor is ignored (the amplitude still vanishes at the axis).
    subroutine test_axis_anchor_pinning()
        integer :: i_anchor
        real(dp) :: arr_in(ns), arr_out(nrho), arr_ref(nrho), anchor
        i_anchor = axis_anchor_index(ns)

        ! m=0: the bare polyfit yields some axis value; pin it to a different one.
        call fill_zernike(0, arr_in)
        call s_to_rho_polyfit(0, ns, nrho, i_anchor, ndeg, arr_in, arr_ref)
        anchor = arr_ref(1) + 0.3_dp
        call s_to_rho_polyfit(0, ns, nrho, i_anchor, ndeg, arr_in, arr_out, anchor=anchor)
        if (abs(arr_out(1) - anchor) > 1.0e-12_dp) then
            print *, 'arr_out(1) =', arr_out(1), ' anchor =', anchor
            error stop 'm=0 anchor not pinned at rho=0'
        end if
        ! The pin must actually move the axis value, not silently no-op.
        if (abs(arr_out(1) - arr_ref(1)) < 1.0e-3_dp) then
            error stop 'anchor had no effect on the axis value'
        end if

        ! m/=0: the anchor must be ignored; the amplitude still vanishes at rho=0.
        call fill_zernike(2, arr_in)
        call s_to_rho_polyfit(2, ns, nrho, i_anchor, ndeg, arr_in, arr_out, anchor=anchor)
        if (abs(arr_out(1)) > 1.0e-13_dp) then
            print *, 'm=2 arr_out(1) with anchor =', arr_out(1)
            error stop 'anchor wrongly applied to m/=0 harmonic'
        end if
    end subroutine test_axis_anchor_pinning

    subroutine fill_zernike(m, arr_in)
        integer, intent(in) :: m
        real(dp), intent(out) :: arr_in(ns)
        integer :: is
        real(dp) :: s
        do is = 1, ns
            s = real(is - 1, dp)/real(ns - 1, dp)
            arr_in(is) = amp*sqrt(s)**m*(1.0_dp + 0.5_dp*s + 0.2_dp*s*s)
        end do
    end subroutine fill_zernike

end program test_axis_polyfit_regularization
