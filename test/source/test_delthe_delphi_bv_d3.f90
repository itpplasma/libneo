!> Validate delthe_delphi_BV_d3 (Boozer angle-map third derivatives) against
!> central finite differences of delthe_delphi_BV_d2 on the LandremanPaul2021 QA
!> reference. The third derivatives feed the analytic d2g_B metric pullback.
program test_delthe_delphi_bv_d3
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
    use boozer_sub, only: get_boozer_coordinates, delthe_delphi_BV_d2, &
                          delthe_delphi_BV_d3
    implicit none

    real(dp), parameter :: reltol = 1.0e-3_dp, floor = 1.0e-7_dp
    character(len=*), parameter :: wout_file = "wout.nc"
    integer, parameter :: n_s = 3, n_ang = 3
    real(dp) :: stor(n_s), vt(n_ang), vp(n_ang)
    integer :: is, ia
    logical :: failed

    use_B_r = .true.
    use_del_tp_B = .true.
    call get_boozer_coordinates(wout_file, radial_spline_order=5, &
                                angular_spline_order=5, grid_refinment=3)

    stor = [0.3_dp, 0.5_dp, 0.7_dp]
    vt = [0.4_dp, 1.7_dp, 2.9_dp]
    vp = [0.2_dp, 0.9_dp, 1.6_dp]

    failed = .false.
    do is = 1, n_s
        do ia = 1, n_ang
            call check_point(stor(is), vt(ia), vp(ia), failed)
        end do
    end do

    if (failed) error stop "delthe_delphi_BV_d3 disagrees with finite differences of d2"
    print *, "delthe_delphi_BV_d3 matches finite differences of delthe_delphi_BV_d2"

contains

    !> Second derivatives (ss,st,sp,tt,tp,pp) of (deltheta, delphi) at (s,t,p).
    subroutine eval_d2(s, t, p, d2)
        real(dp), intent(in) :: s, t, p
        real(dp), intent(out) :: d2(6, 2)
        real(dp) :: dth, dph, g1(3), g2(3), h1(6), h2(6)
        call delthe_delphi_BV_d2(s, t, p, dth, dph, g1, g2, h1, h2)
        d2(:, 1) = h1
        d2(:, 2) = h2
    end subroutine eval_d2

    subroutine d2deriv(s, t, p, h, dir, out)
        real(dp), intent(in) :: s, t, p, h
        integer, intent(in) :: dir
        real(dp), intent(out) :: out(6, 2)
        real(dp) :: dp2(6, 2), dm2(6, 2)
        if (dir == 1) then
            call eval_d2(s + h, t, p, dp2); call eval_d2(s - h, t, p, dm2)
        else if (dir == 2) then
            call eval_d2(s, t + h, p, dp2); call eval_d2(s, t - h, p, dm2)
        else
            call eval_d2(s, t, p + h, dp2); call eval_d2(s, t, p - h, dm2)
        end if
        out = (dp2 - dm2)/(2.0_dp*h)
    end subroutine d2deriv

    subroutine check_point(s, t, p, failed)
        real(dp), intent(in) :: s, t, p
        logical, intent(inout) :: failed
        real(dp) :: dth, dph, g1(3), g2(3), h1(6), h2(6), c1(10), c2(10)
        real(dp) :: d2s(6, 2), d2t(6, 2), d2p(6, 2), fd(10, 2)
        real(dp), parameter :: hs = 2.0e-4_dp, ha = 1.0e-3_dp
        integer :: q

        call delthe_delphi_BV_d3(s, t, p, dth, dph, g1, g2, h1, h2, c1, c2)
        call d2deriv(s, t, p, hs, 1, d2s)
        call d2deriv(s, t, p, ha, 2, d2t)
        call d2deriv(s, t, p, ha, 3, d2p)

        do q = 1, 2
            fd(1, q) = d2s(1, q)   ! sss = d(ss)/ds
            fd(2, q) = d2t(1, q)   ! sst = d(ss)/dt
            fd(3, q) = d2p(1, q)   ! ssp = d(ss)/dp
            fd(4, q) = d2s(4, q)   ! stt = d(tt)/ds
            fd(5, q) = d2s(5, q)   ! stp = d(tp)/ds
            fd(6, q) = d2s(6, q)   ! spp = d(pp)/ds
            fd(7, q) = d2t(4, q)   ! ttt = d(tt)/dt
            fd(8, q) = d2p(4, q)   ! ttp = d(tt)/dp
            fd(9, q) = d2t(6, q)   ! tpp = d(pp)/dt
            fd(10, q) = d2p(6, q)  ! ppp = d(pp)/dp
        end do

        call compare("deltheta", c1, fd(:, 1), failed)
        call compare("delphi", c2, fd(:, 2), failed)
    end subroutine check_point

    subroutine compare(name, analytic, fd, failed)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: analytic(10), fd(10)
        logical, intent(inout) :: failed
        integer :: m
        real(dp) :: denom
        do m = 1, 10
            denom = max(abs(analytic(m)), abs(fd(m)), floor)
            if (abs(analytic(m) - fd(m)) > reltol*denom) then
                print '(a,a,1x,i0,1x,a,es15.7,1x,a,es15.7,1x,es11.3)', name, " d3", &
                    m, "analytic=", analytic(m), "fd=", fd(m), abs(analytic(m) - fd(m))/denom
                failed = .true.
            end if
        end do
    end subroutine compare

end program test_delthe_delphi_bv_d3
