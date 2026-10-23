!> Boozer angle potential along a VMEC theta line. With
!>   vartheta_B = theta + lambda + iota*nu,  varphi_B = varphi + nu
!> the Boozer covariant components are I and G exactly when
!>   (iota*I + G) * d(nu)/d(theta) = B_theta - I - I*d(lambda)/d(theta),
!> where B_theta is the VMEC covariant component. Integrating at fixed varphi,
!>   nu(theta) - nu(0) = [int_0^theta (B_theta - I) dtheta'
!>                        - I*(lambda(theta) - lambda(0))] / (iota*I + G).
!> The right-hand side is evaluated by Gauss-Legendre quadrature of the VMEC
!> field and compared with libneo's nu = delphi_BV. NCSX li383 carries a net
!> toroidal current (I/G ~ 5e-3), so the -I*lambda term is resolved.
program test_boozer_covariant_identity
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_coordinates_mod, only: use_B_r, use_del_tp_B
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord, &
                          delthe_delphi_BV, BOOZER_LEGACY
    use spline_vmec_sub, only: vmec_field
    use new_vmec_stuff_mod, only: nper
    implicit none

    real(dp), parameter :: tol = 1.0e-5_dp
    real(dp), parameter :: s_list(3) = [0.25_dp, 0.5_dp, 0.75_dp]
    real(dp), parameter :: xg(4) = [-0.8611363115940526_dp, -0.3399810435848563_dp, &
                                    0.3399810435848563_dp, 0.8611363115940526_dp]
    real(dp), parameter :: wg(4) = [0.3478548451374538_dp, 0.6521451548625461_dp, &
                                    0.6521451548625461_dp, 0.3478548451374538_dp]
    integer, parameter :: nsub = 400, ncheck = 8, nline = 3
    real(dp) :: pi, s, cur_I, cur_G
    real(dp) :: Ath, Aph, dAth, dAph, d2Aph, d3Aph, dBt, d2Bt, dBp, d2Bp
    real(dp) :: Bm, dBm(3), d2Bm(6), Br, dBr(3), d2Br(6)
    real(dp) :: err
    integer :: k, j

    pi = acos(-1.0_dp)
    use_B_r = .false.
    use_del_tp_B = .false.
    ! This regression integrates the legacy spline field as its oracle.
    call get_boozer_coordinates("wout_li383.nc", construction=BOOZER_LEGACY)

    err = 0.0_dp
    do k = 1, size(s_list)
        s = s_list(k)
        call splint_boozer_coord(s, 0.0_dp, 0.0_dp, 0, Ath, Aph, dAth, dAph, &
                                 d2Aph, d3Aph, cur_I, dBt, d2Bt, cur_G, dBp, &
                                 d2Bp, Bm, dBm, d2Bm, Br, dBr, d2Br)
        call check_line(s, 0.0_dp, 0.0_dp, 1)
        do j = 1, nline
            call check_line(s, 2*pi*(j - 0.4_dp)/nline, 0.0_dp, 2)
        end do
    end do

    print '(A,ES10.3)', 'max |nu(theta) - nu(0) - potential| = ', err
    if (err > tol) then
        print *, 'FAIL: Boozer angle potential misses the -I*lambda term'
        error stop 1
    end if
    print *, 'PASS'

contains

    !> Line from (theta, varphi) along theta (dir=1) or varphi (dir=2) over one
    !> period. The line mean of the covariant component is removed, as libneo
    !> does, so that the check is insensitive to a radial-current residual.
    subroutine check_line(s, th_start, ph_start, dir)
        real(dp), intent(in) :: s, th_start, ph_start
        integer, intent(in) :: dir
        real(dp) :: period, x0, x1, xq, aiota, lam0, lam, b, w(0:nsub), total
        real(dp) :: nu0, nu1, pred
        integer :: i, q

        period = 2*pi
        if (dir == 2) period = 2*pi/nper
        w(0) = 0.0_dp
        do i = 0, nsub - 1
            x0 = period*i/nsub
            x1 = period*(i + 1)/nsub
            w(i + 1) = w(i)
            do q = 1, 4
                xq = 0.5_dp*(x0 + x1) + 0.5_dp*(x1 - x0)*xg(q)
                call field_on_line(s, th_start, ph_start, dir, xq, aiota, lam, b)
                w(i + 1) = w(i + 1) + 0.5_dp*(x1 - x0)*wg(q)*b
            end do
        end do
        total = w(nsub)
        call field_on_line(s, th_start, ph_start, dir, 0.0_dp, aiota, lam0, b)
        call nu_on_line(s, th_start, ph_start, dir, 0.0_dp, nu0)
        do i = nsub/ncheck, nsub, nsub/ncheck
            x1 = period*i/nsub
            call field_on_line(s, th_start, ph_start, dir, x1, aiota, lam, b)
            call nu_on_line(s, th_start, ph_start, dir, x1, nu1)
            pred = (w(i) - total*i/nsub - cur_I*(lam - lam0))/(aiota*cur_I + cur_G)
            err = max(err, abs(nu1 - nu0 - pred))
        end do
    end subroutine check_line

    subroutine field_on_line(s, th_start, ph_start, dir, x, aiota, lam, b)
        real(dp), intent(in) :: s, th_start, ph_start, x
        integer, intent(in) :: dir
        real(dp), intent(out) :: aiota, lam, b
        real(dp) :: At, Ap, dAt, dAp, sqg, dl_ds, dl_dt, dl_dp
        real(dp) :: Bctr_t, Bctr_p, Bcov_s, Bcov_t, Bcov_p, th, ph

        th = th_start
        ph = ph_start
        if (dir == 1) th = th + x
        if (dir == 2) ph = ph + x
        call vmec_field(s, th, ph, At, Ap, dAt, dAp, aiota, sqg, lam, &
                        dl_ds, dl_dt, dl_dp, Bctr_t, Bctr_p, Bcov_s, Bcov_t, Bcov_p)
        if (dir == 1) then
            b = Bcov_t*(1.0_dp + dl_dt)
        else
            b = Bcov_p + Bcov_t*dl_dp
        end if
    end subroutine field_on_line

    subroutine nu_on_line(s, th_start, ph_start, dir, x, nu)
        real(dp), intent(in) :: s, th_start, ph_start, x
        integer, intent(in) :: dir
        real(dp), intent(out) :: nu
        real(dp) :: th, ph, dtheta, ddt(2), ddp(2)

        th = th_start
        ph = ph_start
        if (dir == 1) th = th + x
        if (dir == 2) ph = ph + x
        call delthe_delphi_BV(0, s, th, ph, dtheta, nu, ddt, ddp)
    end subroutine nu_on_line

end program test_boozer_covariant_identity
