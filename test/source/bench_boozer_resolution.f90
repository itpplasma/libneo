!> Self-convergence study of the libneo VMEC->Boozer converter resolution.
!>
!> Builds the Boozer field at several (spline order, multharm) settings and
!> evaluates |B| and sqrt(g) at a fixed set of off-grid points, comparing each
!> setting against a high-resolution reference. Reports the interpolation floor
!> and the grid cost (n_theta*n_phi) so a sensible default multharm can be
!> chosen. Not a pass/fail test; it prints a table.
program bench_boozer_resolution
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_coordinates_mod, only: use_B_r, n_theta_B, n_phi_B
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord, &
                          reset_boozer_batch_splines
    implicit none

    character(len=*), parameter :: wout = "wout.nc"
    integer, parameter :: npts = 200
    real(dp) :: s(npts), th(npts), ph(npts)
    real(dp) :: bmod_ref(npts), sqrtg_ref(npts)
    integer :: j

    ! reference (order 5, fine grid) then the configs to characterise
    integer, parameter :: n_cfg = 7
    integer :: cfg_order(n_cfg), cfg_mh(n_cfg)
    cfg_order = [5, 5, 5, 5, 5, 3, 4]
    cfg_mh    = [3, 4, 5, 6, 7, 3, 3]

    use_B_r = .true.

    do j = 1, npts
        s(j) = 0.05_dp + 0.9_dp*frac(real(j, dp)*0.61803398875_dp)
        th(j) = 2.0_dp*acos(-1.0_dp)*frac(real(j, dp)*0.31830988618_dp)
        ph(j) = 2.0_dp*acos(-1.0_dp)*frac(real(j, dp)*0.13863_dp)
    end do

    ! Reference: order 5, multharm 8.
    call build(5, 8)
    do j = 1, npts
        call eval_bmod_sqrtg(s(j), th(j), ph(j), bmod_ref(j), sqrtg_ref(j))
    end do

    print '(A)', "  order  multharm   n_theta x n_phi    max|dB/B|     rms|dB/B|   max|dsqrtg/sqrtg|"
    print '(A)', "  -------------------------------------------------------------------------------"
    do j = 1, n_cfg
        call report(cfg_order(j), cfg_mh(j))
    end do
    print '(A)', "  (reference: order 5, multharm 8)"

contains

    pure real(dp) function frac(x)
        real(dp), intent(in) :: x
        frac = x - real(int(x), dp)
    end function frac

    subroutine build(order, mh)
        integer, intent(in) :: order, mh
        call reset_boozer_batch_splines()
        call get_boozer_coordinates(wout, radial_spline_order=order, &
                                    angular_spline_order=order, grid_refinment=mh)
    end subroutine build

    subroutine eval_bmod_sqrtg(r, vth, vph, bmod, sqrtg)
        real(dp), intent(in) :: r, vth, vph
        real(dp), intent(out) :: bmod, sqrtg
        real(dp) :: A_theta, A_phi, dA_theta_dr, dA_phi_dr, d2A_phi_dr2, d3A_phi_dr3
        real(dp) :: Bth, dBth, d2Bth, Bph, dBph, d2Bph
        real(dp) :: dBmod(3), d2Bmod(6), Br, dBr(3), d2Br(6)
        real(dp) :: aiota
        call splint_boozer_coord(r, vth, vph, 0, &
                                 A_theta, A_phi, dA_theta_dr, dA_phi_dr, &
                                 d2A_phi_dr2, d3A_phi_dr3, &
                                 Bth, dBth, d2Bth, Bph, dBph, d2Bph, &
                                 bmod, dBmod, d2Bmod, Br, dBr, d2Br)
        aiota = -dA_phi_dr/dA_theta_dr
        ! sqrt(g) up to the constant torflux factor (relative error is unaffected)
        sqrtg = (aiota*Bth + Bph)/bmod**2
    end subroutine eval_bmod_sqrtg

    subroutine report(order, mh)
        integer, intent(in) :: order, mh
        real(dp) :: bmod, sqrtg, db, dg, maxdb, rmsdb, maxdg
        integer :: nt, np_, i
        call build(order, mh)
        nt = n_theta_B
        np_ = n_phi_B
        maxdb = 0.0_dp; rmsdb = 0.0_dp; maxdg = 0.0_dp
        do i = 1, npts
            call eval_bmod_sqrtg(s(i), th(i), ph(i), bmod, sqrtg)
            db = abs(bmod - bmod_ref(i))/abs(bmod_ref(i))
            dg = abs(sqrtg - sqrtg_ref(i))/abs(sqrtg_ref(i))
            maxdb = max(maxdb, db)
            rmsdb = rmsdb + db*db
            maxdg = max(maxdg, dg)
        end do
        rmsdb = sqrt(rmsdb/real(npts, dp))
        print '(I7,I10,I10,A,I4,3X,3ES14.4)', order, mh, nt, " x", np_, maxdb, rmsdb, maxdg
    end subroutine report

end program bench_boozer_resolution
