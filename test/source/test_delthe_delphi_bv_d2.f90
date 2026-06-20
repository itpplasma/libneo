!> Validate the analytic first and second derivatives of the Boozer-side angle
!> map (delthe_delphi_BV_d2) against central finite differences of the routine's
!> own values, on the LandremanPaul2021 QA reference equilibrium. The angle map
!> feeds the curvilinear 6D Boozer metric (boozer_field_metric) dg_B/d2g_B
!> pullback, so its second derivatives must be exact.
program test_delthe_delphi_bv_d2
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
    use boozer_sub, only: get_boozer_coordinates, delthe_delphi_BV_d2
    implicit none

    ! FD truncation on the rho=sqrt(s) composed map dominates the s-direction
    ! second differences; 1e-3 absorbs it while still rejecting any wrong factor.
    real(dp), parameter :: reltol = 1.0e-3_dp, abstol = 1.0e-7_dp
    character(len=*), parameter :: wout_file = "wout.nc"
    integer, parameter :: n_s = 3, n_ang = 3
    real(dp) :: stor(n_s), vt(n_ang), vp(n_ang)
    integer :: is, ia
    logical :: test_failed

    use_B_r = .true.
    use_del_tp_B = .true.
    call get_boozer_coordinates(wout_file, radial_spline_order=5, &
                                angular_spline_order=5, grid_refinment=3)

    stor = [0.3_dp, 0.5_dp, 0.7_dp]
    vt = [0.4_dp, 1.7_dp, 2.9_dp]
    vp = [0.2_dp, 0.9_dp, 1.6_dp]

    test_failed = .false.
    do is = 1, n_s
        do ia = 1, n_ang
            call check_point(stor(is), vt(ia), vp(ia), test_failed)
        end do
    end do

    if (test_failed) then
        error stop "delthe_delphi_BV_d2 derivatives disagree with finite differences"
    end if
    print *, "delthe_delphi_BV_d2 first and second derivatives match finite differences"

contains

    subroutine eval_val(s, t, p, vth, vph)
        real(dp), intent(in) :: s, t, p
        real(dp), intent(out) :: vth, vph
        real(dp) :: g1(3), g2(3), h1(6), h2(6)
        call delthe_delphi_BV_d2(s, t, p, vth, vph, g1, g2, h1, h2)
    end subroutine eval_val

    subroutine eval_shift(q0, k, h, vth, vph)
        real(dp), intent(in) :: q0(3), h
        integer, intent(in) :: k
        real(dp), intent(out) :: vth, vph
        real(dp) :: q(3)
        q = q0
        q(k) = q(k) + h
        call eval_val(q(1), q(2), q(3), vth, vph)
    end subroutine eval_shift

    subroutine eval_shift2(q0, k, l, hk, hl, vth, vph)
        real(dp), intent(in) :: q0(3), hk, hl
        integer, intent(in) :: k, l
        real(dp), intent(out) :: vth, vph
        real(dp) :: q(3)
        q = q0
        q(k) = q(k) + hk
        q(l) = q(l) + hl
        call eval_val(q(1), q(2), q(3), vth, vph)
    end subroutine eval_shift2

    subroutine check_point(s, t, p, test_failed)
        real(dp), intent(in) :: s, t, p
        logical, intent(inout) :: test_failed

        integer, parameter :: idx_i(6) = [1, 1, 1, 2, 2, 3]
        integer, parameter :: idx_j(6) = [1, 2, 3, 2, 3, 3]
        real(dp) :: dth0, dph0, g1(3), g2(3), h1(6), h2(6)
        real(dp) :: hstep(3), q0(3), vth0, vph0
        real(dp) :: fd1_th(3), fd1_ph(3), fd2_th(6), fd2_ph(6)
        real(dp) :: vth_p, vph_p, vth_m, vph_m
        real(dp) :: vth_pp, vph_pp, vth_pm, vph_pm, vth_mp, vph_mp, vth_mm, vph_mm
        integer :: k, l, m

        hstep = [1.0e-3_dp, 1.0e-3_dp, 1.0e-3_dp]
        q0 = [s, t, p]

        call delthe_delphi_BV_d2(s, t, p, dth0, dph0, g1, g2, h1, h2)
        vth0 = dth0
        vph0 = dph0

        do k = 1, 3
            call eval_shift(q0, k, hstep(k), vth_p, vph_p)
            call eval_shift(q0, k, -hstep(k), vth_m, vph_m)
            fd1_th(k) = (vth_p - vth_m)/(2.0_dp*hstep(k))
            fd1_ph(k) = (vph_p - vph_m)/(2.0_dp*hstep(k))
        end do

        do m = 1, 6
            k = idx_i(m)
            l = idx_j(m)
            if (k == l) then
                call eval_shift(q0, k, hstep(k), vth_p, vph_p)
                call eval_shift(q0, k, -hstep(k), vth_m, vph_m)
                fd2_th(m) = (vth_p - 2.0_dp*vth0 + vth_m)/(hstep(k)**2)
                fd2_ph(m) = (vph_p - 2.0_dp*vph0 + vph_m)/(hstep(k)**2)
            else
                call eval_shift2(q0, k, l, hstep(k), hstep(l), vth_pp, vph_pp)
                call eval_shift2(q0, k, l, hstep(k), -hstep(l), vth_pm, vph_pm)
                call eval_shift2(q0, k, l, -hstep(k), hstep(l), vth_mp, vph_mp)
                call eval_shift2(q0, k, l, -hstep(k), -hstep(l), vth_mm, vph_mm)
                fd2_th(m) = (vth_pp - vth_pm - vth_mp + vth_mm) &
                            /(4.0_dp*hstep(k)*hstep(l))
                fd2_ph(m) = (vph_pp - vph_pm - vph_mp + vph_mm) &
                            /(4.0_dp*hstep(k)*hstep(l))
            end if
        end do

        do k = 1, 3
            call compare("dtheta d1", k, g1(k), fd1_th(k), test_failed)
            call compare("dphi d1", k, g2(k), fd1_ph(k), test_failed)
        end do
        do m = 1, 6
            call compare("dtheta d2", m, h1(m), fd2_th(m), test_failed)
            call compare("dphi d2", m, h2(m), fd2_ph(m), test_failed)
        end do
    end subroutine check_point

    subroutine compare(name, idx, analytic, fd, test_failed)
        character(len=*), intent(in) :: name
        integer, intent(in) :: idx
        real(dp), intent(in) :: analytic, fd
        logical, intent(inout) :: test_failed
        real(dp) :: denom

        denom = max(abs(analytic), abs(fd), abstol)
        if (abs(analytic - fd) > reltol*denom) then
            print '(a,1x,i0,1x,a,es16.8,1x,a,es16.8,1x,a,es12.4)', &
                name, idx, "analytic=", analytic, "fd=", fd, &
                "reldiff=", abs(analytic - fd)/denom
            test_failed = .true.
        end if
    end subroutine compare

end program test_delthe_delphi_bv_d2
