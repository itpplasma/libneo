program test_geoflux_flux_table
    ! Accuracy of the geoflux radial table s_tor(psi_pol) against closed forms.
    !
    ! The test writes its own analytic GEQDSK of concentric circular surfaces
    ! about (R0, 0) with constant F = R0*B0 and
    !
    !     psi(rho) = B0/(2c) log(1 + c rho^2/Q0),  c = (QA - Q0)/a^2,
    !
    ! so that psi'(rho) = B0 rho/qm, qm = Q0 + c rho^2.  Every reference value
    ! below follows from these formulas alone, never from a GEQDSK reader:
    !
    !     q(rho)         = qm/sqrt(1 - (rho/R0)^2)
    !     Psi_tor(rho)   = F (R0 - sqrt(R0^2 - rho^2))     (per radian)
    !     s(rho)         = Psi_tor(rho)/Psi_tor(a)
    !
    ! A second, shaped check uses the Solov'ev flux
    !
    !     psi = cs ((R^2 - R0^2)^2/(4 R0^2) + R^2 Z^2/(1.3 R0)^2),  cs = 1.3 B0/(2 Q0)
    !
    ! with reference values from adaptive quadrature (scipy quad, rtol 1e-12)
    ! of F dR dZ/R over the analytic surfaces, not from any GEQDSK reader.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use geoflux_coordinates, only: init_geoflux_coordinates, geoflux_to_cyl, &
        geoflux_get_flux_profiles

    implicit none

    real(dp), parameter :: r0 = 1.65_dp, a = 0.5_dp, b0 = 2.0_dp
    real(dp), parameter :: q0 = 1.2_dp, qa = 3.5_dp
    real(dp), parameter :: c = (qa - q0)/a**2
    ! dq/ds needs the second psi derivative of the 65x65 psirz spline, whose
    ! interpolation error alone is ~2e-4 here (it drops to ~4e-5 on 129x129).
    real(dp), parameter :: tol_s = 1.0e-5_dp, tol_der = 1.0e-4_dp
    real(dp), parameter :: tol_dq = 5.0e-4_dp
    integer, parameter :: nr = 65, nz = 65, nsample = 37
    character(len=*), parameter :: gfile = 'geoflux_flux_table_circ.geqdsk'
    character(len=*), parameter :: sfile = 'geoflux_flux_table_solovev.geqdsk'
    real(dp), parameter :: cs = 1.3_dp*b0/(2.0_dp*q0)
    ! Solov'ev: Psi_tor,edge/psi_b and s at psi_pol/psi_b = 0.35.
    real(dp), parameter :: sol_edge_ref = 1.6142646141116899_dp
    real(dp), parameter :: sol_s_ref = 0.28363823994992976_dp

    real(dp) :: s, rho, q, dq_ds, psi, dpsi_ds, psi_tor_edge, psi_b
    real(dp) :: xgeo(3), xcyl(3), jac(3, 3), d, rho2, drho2_ds, x, qm
    real(dp) :: q_ref, dq_ref, psi_ref, dpsi_ref, drho_ref, s_ref, edge_ref
    real(dp) :: err_s, err_drho, err_q, err_dq, err_psi, err_dpsi, err_edge
    real(dp) :: psi_b_cgs, length_cgs, err_sol_edge, err_sol_psi
    integer :: i, orientation

    d = r0 - sqrt(r0**2 - a**2)
    psi_b = psi_of_rho(a)
    length_cgs = 100.0_dp
    xgeo = 0.0_dp

    edge_ref = b0*r0*d/psi_b
    err_edge = 0.0_dp

    err_s = 0.0_dp
    err_drho = 0.0_dp
    err_q = 0.0_dp
    err_dq = 0.0_dp
    err_psi = 0.0_dp
    err_dpsi = 0.0_dp
    ! Both native poloidal-flux orientations describe the same surfaces.
    do orientation = 1, 2
        call write_geqdsk(gfile, .false., orientation == 2)
        call init_geoflux_coordinates(gfile)
        call geoflux_get_flux_profiles(1.0_dp, q, dq_ds, psi_b_cgs, dpsi_ds, &
            psi_tor_edge)
        if (orientation == 2) then
            if (psi_b_cgs >= 0.0_dp) error stop 'native psi sign changed'
            if (psi_tor_edge >= 0.0_dp) error stop 'native toroidal flux sign changed'
        end if
        err_edge = max(err_edge, abs(abs(psi_tor_edge/psi_b_cgs) - edge_ref)/edge_ref)
        do i = 1, nsample
            s = 0.05_dp + 0.9_dp*real(i - 1, dp)/real(nsample - 1, dp)
            rho2 = r0**2 - (r0 - s*d)**2
            rho = sqrt(rho2)
            drho2_ds = 2.0_dp*d*(r0 - s*d)
            x = rho2/r0**2
            qm = q0 + c*rho2
            q_ref = qm/sqrt(1.0_dp - x)
            dq_ref = (c/sqrt(1.0_dp - x) + 0.5_dp*qm/(r0**2*(1.0_dp - x)**1.5_dp)) &
                *drho2_ds
            psi_ref = psi_of_rho(rho)/psi_b
            dpsi_ref = b0/(2.0_dp*qm)*drho2_ds/psi_b
            drho_ref = 0.5_dp*drho2_ds/rho

            call geoflux_get_flux_profiles(s, q, dq_ds, psi, dpsi_ds, psi_tor_edge)
            err_q = max(err_q, abs(q - q_ref)/q_ref)
            err_dq = max(err_dq, abs(dq_ds - dq_ref)/abs(dq_ref))
            err_psi = max(err_psi, abs(psi/psi_b_cgs - psi_ref))
            err_dpsi = max(err_dpsi, abs(dpsi_ds/psi_b_cgs - dpsi_ref)/dpsi_ref)

            ! Outboard midplane: R - R0 = rho(s) exactly, so the label of the
            ! returned point measures the error of s itself.
            xgeo(1) = s
            call geoflux_to_cyl(xgeo, xcyl, jac)
            rho = xcyl(1)/length_cgs - r0
            s_ref = (r0 - sqrt(r0**2 - rho**2))/d
            err_s = max(err_s, abs(s_ref - s))
            err_drho = max(err_drho, abs(jac(1, 1)/length_cgs - drho_ref)/drho_ref)
        end do
    end do

    ! The Solov'ev grid is taller than the circular one so that every flux
    ! surface lies inside psirz.
    call write_geqdsk(sfile, .true., .false.)
    call init_geoflux_coordinates(sfile)
    call geoflux_get_flux_profiles(1.0_dp, q, dq_ds, psi_b_cgs, dpsi_ds, &
        psi_tor_edge)
    err_sol_edge = abs(abs(psi_tor_edge/psi_b_cgs) - sol_edge_ref)/sol_edge_ref
    call geoflux_get_flux_profiles(sol_s_ref, q, dq_ds, psi, dpsi_ds, psi_tor_edge)
    err_sol_psi = abs(psi/psi_b_cgs - 0.35_dp)

    write (*, '(a,es10.3)') 'Psi_tor edge rel. error  ', err_edge
    write (*, '(a,es10.3)') 's(R) abs. error          ', err_s
    write (*, '(a,es10.3)') 'dR/ds rel. error         ', err_drho
    write (*, '(a,es10.3)') 'psi_pol(s) abs. error    ', err_psi
    write (*, '(a,es10.3)') 'dpsi_pol/ds rel. error   ', err_dpsi
    write (*, '(a,es10.3)') 'q(s) rel. error          ', err_q
    write (*, '(a,es10.3)') 'dq/ds rel. error         ', err_dq

    write (*, '(a,es10.3)') 'Solovev edge flux error  ', err_sol_edge
    write (*, '(a,es10.3)') 'Solovev psi_pol(s) error ', err_sol_psi

    if (err_edge > tol_s) error stop 'toroidal edge flux inaccurate'
    if (err_sol_edge > tol_s) error stop 'Solovev toroidal edge flux inaccurate'
    if (err_sol_psi > tol_s) error stop 'Solovev psi_pol(s) inaccurate'
    if (err_s > tol_s) error stop 's label inaccurate'
    if (err_psi > tol_s) error stop 'psi_pol(s) inaccurate'
    if (err_q > tol_s) error stop 'q(s) inaccurate'
    if (err_drho > tol_der) error stop 'dR/ds inaccurate'
    if (err_dpsi > tol_der) error stop 'dpsi_pol/ds inaccurate'
    if (err_dq > tol_dq) error stop 'dq/ds inaccurate'
    write (*, '(a)') 'PASS'

contains

    pure function psi_of_rho(r) result(p)
        real(dp), intent(in) :: r
        real(dp) :: p

        p = b0/(2.0_dp*c)*log(1.0_dp + c*r**2/q0)
    end function psi_of_rho

    pure function psi_rz(r, z, solovev) result(p)
        real(dp), intent(in) :: r, z
        logical, intent(in) :: solovev
        real(dp) :: p

        if (solovev) then
            p = cs*((r**2 - r0**2)**2/(4.0_dp*r0**2) + r**2*z**2/(1.3_dp*r0)**2)
        else
            p = psi_of_rho(hypot(r - r0, z))
        end if
    end function psi_rz

    subroutine write_geqdsk(path, solovev, reverse_psi)
        character(len=*), intent(in) :: path
        logical, intent(in) :: solovev, reverse_psi

        integer, parameter :: nb = 129
        real(dp) :: rdim, zdim, rleft, rr(nr), zz(nz), psirz(nr, nz)
        real(dp) :: qpsi(nr), fpol(nr), zeros(nr), th
        real(dp) :: lcfs(2, nb), lim(2, nb), sibry, psi_sign
        integer :: u, i, j

        psi_sign = merge(-1.0_dp, 1.0_dp, reverse_psi)
        rdim = merge(3.4_dp, 3.2_dp, solovev)*a
        zdim = merge(3.6_dp, 3.2_dp, solovev)*a
        rleft = r0 - 0.5_dp*rdim
        do i = 1, nr
            rr(i) = rleft + rdim*real(i - 1, dp)/real(nr - 1, dp)
        end do
        do j = 1, nz
            zz(j) = -0.5_dp*zdim + zdim*real(j - 1, dp)/real(nz - 1, dp)
        end do
        do j = 1, nz
            do i = 1, nr
                psirz(i, j) = psi_sign*psi_rz(rr(i), zz(j), solovev)
            end do
        end do
        sibry = psi_sign*psi_rz(r0 + a, 0.0_dp, solovev)
        ! qpsi is deliberately constant and the boundary record is a plain
        ! circle: geoflux must use neither.
        qpsi = q0
        fpol = r0*b0
        zeros = 0.0_dp
        do i = 1, nb
            th = 2.0_dp*acos(-1.0_dp)*real(i - 1, dp)/real(nb - 1, dp)
            lcfs(1, i) = r0 + a*cos(th)
            lcfs(2, i) = a*sin(th)
            lim(1, i) = r0 + 1.5_dp*a*cos(th)
            lim(2, i) = 1.5_dp*a*sin(th)
        end do

        open (newunit=u, file=path, status='replace', action='write')
        write (u, '(a48,3i4)') 'libneo analytic', 0, nr, nz
        write (u, '(5es16.9)') rdim, zdim, r0, rleft, 0.0_dp
        write (u, '(5es16.9)') r0, 0.0_dp, 0.0_dp, sibry, b0
        write (u, '(5es16.9)') psi_sign*1.0e6_dp, 0.0_dp, 0.0_dp, r0, 0.0_dp
        write (u, '(5es16.9)') 0.0_dp, 0.0_dp, sibry, 0.0_dp, 0.0_dp
        write (u, '(5es16.9)') fpol
        write (u, '(5es16.9)') zeros
        write (u, '(5es16.9)') zeros
        write (u, '(5es16.9)') zeros
        write (u, '(5es16.9)') psirz
        write (u, '(5es16.9)') qpsi
        write (u, '(2i5)') nb, nb
        write (u, '(5es16.9)') lcfs
        write (u, '(5es16.9)') lim
        close (u)
    end subroutine write_geqdsk

end program test_geoflux_flux_table
