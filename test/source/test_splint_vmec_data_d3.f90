!> Validate splint_vmec_data_d3 (third derivatives of the VMEC coordinate map
!> R, Z, lambda) against central finite differences of splint_vmec_data_d2 on
!> the LandremanPaul2021 QA reference equilibrium. The third derivatives feed
!> the analytic second derivative of the curvilinear metric (boozer_field_metric
!> d2g_B), so they must match the validated second-derivative routine.
program test_splint_vmec_data_d3
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use new_vmec_stuff_mod, only: netcdffile
    use spline_vmec_sub, only: spline_vmec_data, splint_vmec_data_d2, &
                               splint_vmec_data_d3
    implicit none

    real(dp), parameter :: reltol = 5.0e-3_dp, floor = 1.0e-3_dp
    integer, parameter :: n_pts = 4
    real(dp) :: stor(n_pts), tt(n_pts), pp(n_pts)
    integer :: ip
    logical :: failed

    netcdffile = "wout.nc"
    call spline_vmec_data

    stor = [0.3_dp, 0.5_dp, 0.6_dp, 0.7_dp]
    tt = [0.4_dp, 1.7_dp, 2.9_dp, 0.9_dp]
    pp = [0.2_dp, 0.9_dp, 1.6_dp, 2.4_dp]

    failed = .false.
    do ip = 1, n_pts
        call check_point(stor(ip), tt(ip), pp(ip), failed)
    end do

    if (failed) error stop "splint_vmec_data_d3 disagrees with finite differences of d2"
    print *, "splint_vmec_data_d3 matches finite differences of splint_vmec_data_d2"

contains

    !> Second derivatives (ss, st, sp, tt, tp, pp) of R, Z, lambda, stacked as
    !> d2(field, comp) with field 1=R, 2=Z, 3=lambda.
    subroutine eval_d2(s, t, p, d2)
        real(dp), intent(in) :: s, t, p
        real(dp), intent(out) :: d2(3, 6)
        real(dp) :: A_phi, A_theta, dA_phi_ds, dA_theta_ds, aiota, R, Z, alam
        real(dp) :: dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, dl_ds, dl_dt, dl_dp
        real(dp) :: d2R(6), d2Z(6), d2l(6)
        call splint_vmec_data_d2(s, t, p, A_phi, A_theta, dA_phi_ds, dA_theta_ds, &
            aiota, R, Z, alam, dR_ds, dR_dt, dR_dp, dZ_ds, dZ_dt, dZ_dp, &
            dl_ds, dl_dt, dl_dp, d2R, d2Z, d2l)
        d2(1, :) = d2R
        d2(2, :) = d2Z
        d2(3, :) = d2l
    end subroutine eval_d2

    subroutine check_point(s, t, p, failed)
        real(dp), intent(in) :: s, t, p
        logical, intent(inout) :: failed
        real(dp) :: hs, ht, hp
        real(dp) :: d2sp(3, 6), d2tp(3, 6), d2pp(3, 6)
        real(dp) :: dp2(3, 6), dm2(3, 6)
        real(dp) :: d3R(10), d3Z(10), d3l(10), fd(3, 10)
        integer :: f

        hs = 2.0e-4_dp
        ht = 1.0e-3_dp
        hp = 1.0e-3_dp

        ! d(d2)/ds, d(d2)/dt, d(d2)/dp by central differences.
        call eval_d2(s + hs, t, p, dp2); call eval_d2(s - hs, t, p, dm2)
        d2sp = (dp2 - dm2)/(2.0_dp*hs)
        call eval_d2(s, t + ht, p, dp2); call eval_d2(s, t - ht, p, dm2)
        d2tp = (dp2 - dm2)/(2.0_dp*ht)
        call eval_d2(s, t, p + hp, dp2); call eval_d2(s, t, p - hp, dm2)
        d2pp = (dp2 - dm2)/(2.0_dp*hp)

        ! Assemble the 10 third derivatives (sss,sst,ssp,stt,stp,spp,ttt,ttp,tpp,ppp)
        ! from d(d2_ab)/dc. d2 packed (ss=1, st=2, sp=3, tt=4, tp=5, pp=6).
        do f = 1, 3
            fd(f, 1) = d2sp(f, 1)   ! sss = d(ss)/ds
            fd(f, 2) = d2tp(f, 1)   ! sst = d(ss)/dt
            fd(f, 3) = d2pp(f, 1)   ! ssp = d(ss)/dp
            fd(f, 4) = d2sp(f, 4)   ! stt = d(tt)/ds
            fd(f, 5) = d2sp(f, 5)   ! stp = d(tp)/ds
            fd(f, 6) = d2sp(f, 6)   ! spp = d(pp)/ds
            fd(f, 7) = d2tp(f, 4)   ! ttt = d(tt)/dt
            fd(f, 8) = d2pp(f, 4)   ! ttp = d(tt)/dp
            fd(f, 9) = d2tp(f, 6)   ! tpp = d(pp)/dt
            fd(f, 10) = d2pp(f, 6)  ! ppp = d(pp)/dp
        end do

        call splint_vmec_data_d3(s, t, p, d3R, d3Z, d3l)
        call compare_field("R", d3R, fd(1, :), failed)
        call compare_field("Z", d3Z, fd(2, :), failed)
        call compare_field("lambda", d3l, fd(3, :), failed)
    end subroutine check_point

    subroutine compare_field(name, analytic, fd, failed)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: analytic(10), fd(10)
        logical, intent(inout) :: failed
        integer :: m
        real(dp) :: denom

        do m = 1, 10
            denom = max(abs(analytic(m)), abs(fd(m)), floor)
            if (abs(analytic(m) - fd(m)) > reltol*denom) then
                print '(a,a,1x,i0,1x,a,es15.7,1x,a,es15.7,1x,a,es11.3)', &
                    name, " d3", m, "analytic=", analytic(m), "fd=", fd(m), &
                    "reldiff=", abs(analytic(m) - fd(m))/denom
                failed = .true.
            end if
        end do
    end subroutine compare_field

end program test_splint_vmec_data_d3
