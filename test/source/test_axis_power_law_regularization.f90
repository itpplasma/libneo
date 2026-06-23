program test_axis_power_law_regularization
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use new_vmec_stuff_mod, only: ns_s, s_axis_heal, rho_axis_heal
    use spline_vmec_sub, only: s_to_rho_power_law, axis_anchor_index
    implicit none

    integer, parameter :: ns = 100, nrho = 100
    real(dp), parameter :: amp = 0.75_dp

    ns_s = 5

    call test_anchor_index()
    call test_axis_vanishing()
    call test_power_law_reproduction()
    call test_anchor_continuity()

    print *, 'axis power-law regularization: all checks passed'

contains

    !> The innermost reliable surface index used by the rho**|m| continuation.
    subroutine test_anchor_index()
        s_axis_heal = 1.0e-2_dp
        rho_axis_heal = -1.0_dp
        if (axis_anchor_index(ns) /= 2) then
            print *, 'i_anchor(s_axis_heal=0.01) =', axis_anchor_index(ns), 'expected 2'
            error stop 'axis_anchor_index mismatch for s_axis_heal=0.01'
        end if
        s_axis_heal = 9.0e-2_dp
        if (axis_anchor_index(ns) /= 10) then
            print *, 'i_anchor(s_axis_heal=0.09) =', axis_anchor_index(ns), 'expected 10'
            error stop 'axis_anchor_index mismatch for s_axis_heal=0.09'
        end if
        rho_axis_heal = 0.3_dp
        if (axis_anchor_index(ns) /= 10) then
            print *, 'i_anchor(rho_axis_heal=0.3) =', axis_anchor_index(ns), 'expected 10'
            error stop 'deprecated rho_axis_heal conversion mismatch'
        end if
        s_axis_heal = 1.0e-2_dp
        rho_axis_heal = -1.0_dp
    end subroutine test_anchor_index

    !> A harmonic with m>0 must vanish at the magnetic axis (rho=0).
    subroutine test_axis_vanishing()
        integer :: i_anchor
        real(dp) :: arr_in(ns), arr_out(nrho)

        s_axis_heal = 1.0e-2_dp
        rho_axis_heal = -1.0_dp
        i_anchor = axis_anchor_index(ns)
        call fill_power_law(2, arr_in)
        call s_to_rho_power_law(2, ns, nrho, i_anchor, arr_in, arr_out)
        if (abs(arr_out(1)) > 1.0e-13_dp) then
            print *, 'arr_out at rho=0 =', arr_out(1)
            error stop 'm>0 harmonic does not vanish at the axis'
        end if
    end subroutine test_axis_vanishing

    !> An exact c(rho) = amp*rho**m input is an even-analytic reduced
    !> function (c/rho**m = amp), so the resampled harmonic must reproduce
    !> amp*rho**m on the whole rho grid, including below the anchor.
    subroutine test_power_law_reproduction()
        integer, parameter :: ms(3) = [1, 2, 4]
        integer :: im, m, i_anchor, irho
        real(dp) :: arr_in(ns), arr_out(nrho)
        real(dp) :: rho, expected, err, max_err

        s_axis_heal = 1.0e-2_dp
        rho_axis_heal = -1.0_dp
        i_anchor = axis_anchor_index(ns)
        do im = 1, size(ms)
            m = ms(im)
            call fill_power_law(m, arr_in)
            call s_to_rho_power_law(m, ns, nrho, i_anchor, arr_in, arr_out)
            max_err = 0.0_dp
            do irho = 1, nrho
                rho = real(irho - 1, dp)/real(nrho - 1, dp)
                expected = amp*rho**m
                err = abs(arr_out(irho) - expected)
                max_err = max(max_err, err)
            end do
            if (max_err > 1.0e-10_dp) then
                print *, 'm =', m, 'max |arr_out - amp*rho**m| =', max_err
                error stop 'power-law harmonic not reproduced (see below-anchor branch)'
            end if
        end do
    end subroutine test_power_law_reproduction

    !> The continuation must be continuous across the anchor surface.
    subroutine test_anchor_continuity()
        integer :: i_anchor, irho_below, irho_above, irho
        real(dp) :: arr_in(ns), arr_out(nrho)
        real(dp) :: rho, rho_anchor, jump

        s_axis_heal = 1.0e-2_dp
        rho_axis_heal = -1.0_dp
        i_anchor = axis_anchor_index(ns)
        rho_anchor = sqrt(real(i_anchor - 1, dp)/real(ns - 1, dp))
        call fill_power_law(2, arr_in)
        call s_to_rho_power_law(2, ns, nrho, i_anchor, arr_in, arr_out)

        irho_below = 0
        irho_above = 0
        do irho = 1, nrho
            rho = real(irho - 1, dp)/real(nrho - 1, dp)
            if (rho < rho_anchor) irho_below = irho
            if (rho >= rho_anchor .and. irho_above == 0) irho_above = irho
        end do

        jump = abs(arr_out(irho_above) - arr_out(irho_below))
        if (jump > 5.0e-3_dp) then
            print *, 'rho_anchor =', rho_anchor
            print *, 'arr_out just below anchor =', arr_out(irho_below)
            print *, 'arr_out just above anchor =', arr_out(irho_above)
            print *, 'jump across anchor =', jump
            error stop 'discontinuity in rho**|m| continuation at the anchor'
        end if
    end subroutine test_anchor_continuity

    !> Sample c(s) = amp*s**(m/2) = amp*rho**m on the uniform s grid.
    subroutine fill_power_law(m, arr_in)
        integer, intent(in) :: m
        real(dp), intent(out) :: arr_in(ns)
        integer :: is
        real(dp) :: s

        do is = 1, ns
            s = real(is - 1, dp)/real(ns - 1, dp)
            arr_in(is) = amp*sqrt(s)**m
        end do
    end subroutine fill_power_law

end program test_axis_power_law_regularization
