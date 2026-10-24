!> One flux surface of the VMEC -> Boozer transform as a function of its
!> inputs, with its tangent (JVP) and adjoint (VJP).
!>
!> Inputs are the field data on the equidistant VMEC angle grid
!> (theta_V, varphi_V) of one surface: iota, lambda, the Boozer generating
!> components B_theta(1 + lambda_theta) and B_varphi + B_theta lambda_varphi,
!> and further quantities (slot 3 = |B|) that are carried along to the
!> equidistant Boozer grid. boozer_surface_transform is the primal used by
!> boozer_sub/compute_boozer_data. The linearised routines differentiate the
!> same discrete algorithm:
!>  - the flux-surface averages I and G and the periodic-spline integration of
!>    the generating components are linear; their tangent applies explicit
!>    operator matrices and their adjoint the transposed matrices;
!>  - the re-interpolation to Boozer angles is Lagrange interpolation at
!>    moving nodes x_k = x_V + Delta(x_k). Its exact derivative uses
!>    d p(x)/d x_j = -L_j(x) p'(x_j), so the tangent is the interpolation of
!>    the materially corrected node tangents q_dot_k - x_dot_k p'(x_k); this is
!>    the discrete counterpart of the implicit function theorem for the angle
!>    inversion theta_V + Delta(theta_V) = theta_B.
!> The stencil chosen by the binary search is a discrete choice. It is frozen
!> in the derivative, so the derivative is exact only away from stencil switch
!> points, where the primal itself jumps at the interpolation-error level. For
!> stellarator-symmetric data the angle shifts vanish at the symmetry points,
!> a moving node then lands exactly on its target and the stencil is a tie:
!> the value is continuous there but the two one-sided derivatives differ by
!> the interpolation error of p'(x), and the tangent returns the one-sided
!> derivative of the stencil the binary search picks.
module boozer_surface
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spl_three_to_five_sub, only: spl_per
    use binsrc_sub, only: binsrc
    use plag_coeff_sub, only: plag_coeff

    implicit none
    private

    public :: boozer_angle_grid_t, boozer_surface_lin_t
    public :: boozer_surface_transform, boozer_surface_lin_init
    public :: boozer_surface_jvp, boozer_surface_vjp
    public :: periodic_cumulative_integral

    !> Equidistant periodic angle grid shared by VMEC and Boozer angles. The
    !> last node duplicates the first one (theta = 2 pi, varphi = 2 pi/nper).
    type :: boozer_angle_grid_t
        integer :: n_theta = 0
        integer :: n_phi = 0
        integer :: ns_tp = 5
        real(dp) :: h_theta = 0.0_dp
        real(dp) :: h_phi = 0.0_dp
    end type boozer_angle_grid_t

    !> Linear operators of the surface transform, built once per grid.
    !> k_theta and k_phi map nodal values to the detrended cumulative periodic
    !> spline integral starting at zero.
    type :: boozer_surface_lin_t
        type(boozer_angle_grid_t) :: grid
        real(dp), allocatable :: k_theta(:, :)
        real(dp), allocatable :: k_phi(:, :)
    end type boozer_surface_lin_t

    !> Forward pass shared by the JVP and the VJP; it keeps the stencils,
    !> weights and nodal derivatives of both re-interpolation sweeps.
    type :: surface_tape_t
        real(dp) :: denomjac, gbeg, bcov_t, bcov_p
        real(dp), allocatable :: raw(:, :), dphi(:, :), dtheta(:, :)
        real(dp), allocatable :: q0(:, :, :), q1(:, :, :), q2(:, :, :)
        integer, allocatable :: ib1(:, :), ib2(:, :)
        real(dp), allocatable :: w1(:, :, :), w2(:, :, :)
        real(dp), allocatable :: p1(:, :, :, :), p2(:, :, :, :)
    end type surface_tape_t

contains

    !> Detrended cumulative integral of a periodic spline through f, starting
    !> at the value start. This is the exact operation sequence of the legacy
    !> converter, so the primal stays bit-identical.
    subroutine periodic_cumulative_integral(ns_tp, n, h, wint, f, start, res)
        integer, intent(in) :: ns_tp, n
        real(dp), intent(in) :: h, wint(0:ns_tp), f(n), start
        real(dp), intent(out) :: res(n)
        real(dp) :: splcoe(0:ns_tp, n), aper
        integer :: i

        splcoe(0, :) = f
        call spl_per(ns_tp, n, h, splcoe)
        res(1) = start
        do i = 1, n - 1
            res(i + 1) = res(i) + sum(wint*splcoe(:, i))
        end do
        aper = (res(n) - res(1))/real(n - 1, dp)
        do i = 2, n
            res(i) = res(i) - aper*real(i - 1, dp)
        end do
    end subroutine periodic_cumulative_integral

    subroutine integration_weights(ns_tp, h, wint)
        integer, intent(in) :: ns_tp
        real(dp), intent(in) :: h
        real(dp), intent(out) :: wint(0:ns_tp)
        integer :: i

        do i = 0, ns_tp
            wint(i) = h**(i + 1)/real(i + 1, dp)
        end do
    end subroutine integration_weights

    !> Equidistant nodes extended by one period on both sides.
    subroutine extended_nodes(n, h, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: h
        real(dp), intent(out) :: x(2 - n:2*n - 1)
        real(dp) :: per
        integer :: i

        do i = 1, n
            x(i) = real(i - 1, dp)*h
        end do
        per = real(n - 1, dp)*h
        x(2 - n:0) = x(1:n - 1) - per
        x(n + 1:2*n - 1) = x(2:n) + per
    end subroutine extended_nodes

    !> Primal transform of one surface.
    !>
    !> perqua(3:, :, :) holds the carried quantities on the VMEC grid on entry;
    !> on exit perqua(:, :, :) holds (Delta theta, Delta varphi, carried
    !> quantities) on the Boozer grid. deltheta_V and delphi_V return the angle
    !> differences theta_B - theta and varphi_B - varphi on the VMEC grid.
    subroutine boozer_surface_transform(grid, aiota, alam, bcov_theta_V, &
                                        bcov_varphi_V, perqua, &
                                        bcov_vartheta_B, bcov_varphi_B, &
                                        deltheta_V, delphi_V)
        type(boozer_angle_grid_t), intent(in) :: grid
        real(dp), intent(in) :: aiota
        real(dp), intent(in) :: alam(:, :), bcov_theta_V(:, :), bcov_varphi_V(:, :)
        real(dp), intent(inout) :: perqua(:, :, :)
        real(dp), intent(out) :: bcov_vartheta_B, bcov_varphi_B
        real(dp), intent(out) :: deltheta_V(:, :), delphi_V(:, :)

        integer :: n_theta, n_phi, ns_tp, nqua, npoilag, nshift
        integer :: i, i_theta, i_phi, ibeg, iend
        real(dp) :: gridcellnum, denomjac, Gbeg, start
        real(dp), allocatable :: wint_t(:), wint_p(:), coef(:, :)
        real(dp), allocatable :: theta_V(:), theta_B(:), phi_V(:), phi_B(:)
        real(dp), allocatable :: perqua_t(:, :), perqua_p(:, :)

        n_theta = grid%n_theta
        n_phi = grid%n_phi
        ns_tp = grid%ns_tp
        nqua = size(perqua, 1)
        npoilag = ns_tp + 1
        nshift = npoilag/2
        gridcellnum = real((n_theta - 1)*(n_phi - 1), dp)

        allocate (wint_t(0:ns_tp), wint_p(0:ns_tp), coef(0:0, npoilag))
        allocate (theta_V(2 - n_theta:2*n_theta - 1))
        allocate (theta_B(2 - n_theta:2*n_theta - 1))
        allocate (phi_V(2 - n_phi:2*n_phi - 1), phi_B(2 - n_phi:2*n_phi - 1))
        allocate (perqua_t(nqua, 2 - n_theta:2*n_theta - 1))
        allocate (perqua_p(nqua, 2 - n_phi:2*n_phi - 1))

        call integration_weights(ns_tp, grid%h_theta, wint_t)
        call integration_weights(ns_tp, grid%h_phi, wint_p)
        call extended_nodes(n_theta, grid%h_theta, theta_V)
        call extended_nodes(n_phi, grid%h_phi, phi_V)

        bcov_vartheta_B = sum(bcov_theta_V(2:n_theta, 2:n_phi))/gridcellnum
        bcov_varphi_B = sum(bcov_varphi_V(2:n_theta, 2:n_phi))/gridcellnum

        denomjac = 1.0_dp/(aiota*bcov_vartheta_B + bcov_varphi_B)
        Gbeg = bcov_vartheta_B*denomjac*alam(1, 1)

        call periodic_cumulative_integral(ns_tp, n_theta, grid%h_theta, wint_t, &
                                          bcov_theta_V(:, 1), 0.0_dp, &
                                          delphi_V(:, 1))
        do i_theta = 1, n_theta
            start = delphi_V(i_theta, 1)
            call periodic_cumulative_integral(ns_tp, n_phi, grid%h_phi, wint_p, &
                                              bcov_varphi_V(i_theta, :), start, &
                                              delphi_V(i_theta, :))
        end do

        ! Boozer condition (iota B_vartheta + B_varphi) dG = dw - B_vartheta dlambda
        delphi_V = denomjac*(delphi_V - bcov_vartheta_B*alam) + Gbeg
        deltheta_V = aiota*delphi_V + alam

        do i_phi = 1, n_phi
            perqua_t(1, 1:n_theta) = deltheta_V(:, i_phi)
            perqua_t(2, 1:n_theta) = delphi_V(:, i_phi)
            perqua_t(3:nqua, 1:n_theta) = perqua(3:nqua, :, i_phi)
            perqua_t(:, 2 - n_theta:0) = perqua_t(:, 1:n_theta - 1)
            perqua_t(:, n_theta + 1:2*n_theta - 1) = perqua_t(:, 2:n_theta)
            theta_B = theta_V + perqua_t(1, :)
            do i_theta = 1, n_theta
                call binsrc(theta_B, 2 - n_theta, 2*n_theta - 1, &
                            theta_V(i_theta), i)
                ibeg = i - nshift
                iend = ibeg + ns_tp
                call plag_coeff(npoilag, 0, theta_V(i_theta), &
                                theta_B(ibeg:iend), coef)
                perqua(:, i_theta, i_phi) = matmul(perqua_t(:, ibeg:iend), &
                                                   coef(0, :))
            end do
        end do

        do i_theta = 1, n_theta
            perqua_p(:, 1:n_phi) = perqua(:, i_theta, :)
            perqua_p(:, 2 - n_phi:0) = perqua_p(:, 1:n_phi - 1)
            perqua_p(:, n_phi + 1:2*n_phi - 1) = perqua_p(:, 2:n_phi)
            phi_B = phi_V + perqua_p(2, :)
            do i_phi = 1, n_phi
                call binsrc(phi_B, 2 - n_phi, 2*n_phi - 1, phi_V(i_phi), i)
                ibeg = i - nshift
                iend = ibeg + ns_tp
                call plag_coeff(npoilag, 0, phi_V(i_phi), phi_B(ibeg:iend), coef)
                perqua(:, i_theta, i_phi) = matmul(perqua_p(:, ibeg:iend), &
                                                   coef(0, :))
            end do
        end do
    end subroutine boozer_surface_transform

    !> Build the integration operators by applying the primal integral to
    !> unit vectors; the integral is linear in its data, so these matrices are
    !> its exact Jacobian.
    subroutine boozer_surface_lin_init(grid, lin)
        type(boozer_angle_grid_t), intent(in) :: grid
        type(boozer_surface_lin_t), intent(out) :: lin

        lin%grid = grid
        call integral_matrix(grid%ns_tp, grid%n_theta, grid%h_theta, lin%k_theta)
        call integral_matrix(grid%ns_tp, grid%n_phi, grid%h_phi, lin%k_phi)
    end subroutine boozer_surface_lin_init

    subroutine integral_matrix(ns_tp, n, h, k)
        integer, intent(in) :: ns_tp, n
        real(dp), intent(in) :: h
        real(dp), allocatable, intent(out) :: k(:, :)
        real(dp) :: wint(0:ns_tp), e(n)
        integer :: j

        allocate (k(n, n))
        call integration_weights(ns_tp, h, wint)
        do j = 1, n
            e = 0.0_dp
            e(j) = 1.0_dp
            call periodic_cumulative_integral(ns_tp, n, h, wint, e, 0.0_dp, k(:, j))
        end do
    end subroutine integral_matrix

    !> Base index of an extended periodic index (last node duplicates first).
    pure integer function base_index(e, n)
        integer, intent(in) :: e, n

        base_index = e
        if (e < 1) base_index = e + n - 1
        if (e > n) base_index = e - n + 1
    end function base_index

    !> Lagrange weights at x for the nodes xs and the nodal derivative matrix
    !> d(k, j) = L_j'(xs(k)).
    pure subroutine lagrange_weights(x, xs, w, d)
        real(dp), intent(in) :: x, xs(:)
        real(dp), intent(out) :: w(size(xs)), d(size(xs), size(xs))
        real(dp) :: bw(size(xs))
        integer :: n, j, k, l

        n = size(xs)
        do j = 1, n
            w(j) = 1.0_dp
            bw(j) = 1.0_dp
            do l = 1, n
                if (l == j) cycle
                w(j) = w(j)*(x - xs(l))/(xs(j) - xs(l))
                bw(j) = bw(j)*(xs(j) - xs(l))
            end do
        end do
        do k = 1, n
            d(k, k) = 0.0_dp
            do j = 1, n
                if (j == k) cycle
                d(k, j) = (bw(k)/bw(j))/(xs(k) - xs(j))
                d(k, k) = d(k, k) + 1.0_dp/(xs(k) - xs(j))
            end do
        end do
    end subroutine lagrange_weights

    !> Moving-node interpolation along one periodic direction: values q(:, 1:n)
    !> at nodes x_V + q(idx, :) are interpolated to the targets x_V(1:n).
    !> Returns the stencil start (extended index), weights and nodal
    !> derivatives p'(x_k) of the interpolant for every target.
    subroutine moving_interp(q, idx, xv, ns_tp, res, ibeg, w, pder)
        real(dp), intent(in) :: q(:, :)
        integer, intent(in) :: idx, ns_tp
        real(dp), intent(in) :: xv(2 - size(q, 2):)
        real(dp), intent(out) :: res(size(q, 1), size(q, 2))
        integer, intent(out) :: ibeg(size(q, 2))
        real(dp), intent(out) :: w(ns_tp + 1, size(q, 2))
        real(dp), intent(out) :: pder(size(q, 1), ns_tp + 1, size(q, 2))
        integer :: n, t, k, j, i, e
        real(dp) :: xb(2 - size(q, 2):2*size(q, 2) - 1), d(ns_tp + 1, ns_tp + 1)
        real(dp) :: qs(size(q, 1), ns_tp + 1)

        n = size(q, 2)
        do e = 2 - n, 2*n - 1
            xb(e) = xv(e) + q(idx, base_index(e, n))
        end do
        do t = 1, n
            call binsrc(xb, 2 - n, 2*n - 1, xv(t), i)
            ibeg(t) = i - (ns_tp + 1)/2
            do k = 1, ns_tp + 1
                qs(:, k) = q(:, base_index(ibeg(t) + k - 1, n))
            end do
            call lagrange_weights(xv(t), xb(ibeg(t):ibeg(t) + ns_tp), w(:, t), d)
            res(:, t) = matmul(qs, w(:, t))
            do k = 1, ns_tp + 1
                pder(:, k, t) = 0.0_dp
                do j = 1, ns_tp + 1
                    pder(:, k, t) = pder(:, k, t) + d(k, j)*qs(:, j)
                end do
            end do
        end do
    end subroutine moving_interp

    subroutine moving_interp_jvp(qd, idx, ns_tp, ibeg, w, pder, resd)
        real(dp), intent(in) :: qd(:, :)
        integer, intent(in) :: idx, ns_tp
        integer, intent(in) :: ibeg(:)
        real(dp), intent(in) :: w(:, :), pder(:, :, :)
        real(dp), intent(out) :: resd(size(qd, 1), size(qd, 2))
        integer :: n, t, k, b

        n = size(qd, 2)
        do t = 1, n
            resd(:, t) = 0.0_dp
            do k = 1, ns_tp + 1
                b = base_index(ibeg(t) + k - 1, n)
                resd(:, t) = resd(:, t) + w(k, t)*(qd(:, b) - qd(idx, b)*pder(:, k, t))
            end do
        end do
    end subroutine moving_interp_jvp

    subroutine moving_interp_vjp(resb, idx, ns_tp, ibeg, w, pder, qb)
        real(dp), intent(in) :: resb(:, :)
        integer, intent(in) :: idx, ns_tp
        integer, intent(in) :: ibeg(:)
        real(dp), intent(in) :: w(:, :), pder(:, :, :)
        real(dp), intent(inout) :: qb(:, :)
        integer :: n, t, k, b
        real(dp) :: vb(size(resb, 1))

        n = size(resb, 2)
        do t = 1, n
            do k = 1, ns_tp + 1
                b = base_index(ibeg(t) + k - 1, n)
                vb = w(k, t)*resb(:, t)
                qb(:, b) = qb(:, b) + vb
                qb(idx, b) = qb(idx, b) - dot_product(vb, pder(:, k, t))
            end do
        end do
    end subroutine moving_interp_vjp

    subroutine surface_forward(lin, aiota, alam, bt, bp, carried, tp)
        type(boozer_surface_lin_t), intent(in) :: lin
        real(dp), intent(in) :: aiota, alam(:, :), bt(:, :), bp(:, :)
        real(dp), intent(in) :: carried(:, :, :)
        type(surface_tape_t), intent(out) :: tp
        integer :: nt, np, nq, ns_tp, i, j
        real(dp) :: cells, tv(2 - lin%grid%n_theta:2*lin%grid%n_theta - 1)
        real(dp) :: pv(2 - lin%grid%n_phi:2*lin%grid%n_phi - 1)

        nt = lin%grid%n_theta
        np = lin%grid%n_phi
        ns_tp = lin%grid%ns_tp
        nq = size(carried, 1) + 2
        cells = real((nt - 1)*(np - 1), dp)
        call extended_nodes(nt, lin%grid%h_theta, tv)
        call extended_nodes(np, lin%grid%h_phi, pv)

        tp%bcov_t = sum(bt(2:nt, 2:np))/cells
        tp%bcov_p = sum(bp(2:nt, 2:np))/cells
        tp%denomjac = 1.0_dp/(aiota*tp%bcov_t + tp%bcov_p)
        tp%gbeg = tp%bcov_t*tp%denomjac*alam(1, 1)

        allocate (tp%raw(nt, np))
        tp%raw(:, 1) = matmul(lin%k_theta, bt(:, 1))
        do i = 1, nt
            tp%raw(i, :) = tp%raw(i, 1) + matmul(lin%k_phi, bp(i, :))
        end do
        allocate (tp%dphi(nt, np), tp%dtheta(nt, np))
        tp%dphi = tp%denomjac*(tp%raw - tp%bcov_t*alam) + tp%gbeg
        tp%dtheta = aiota*tp%dphi + alam

        allocate (tp%q0(nq, nt, np), tp%q1(nq, nt, np), tp%q2(nq, nt, np))
        tp%q0(1, :, :) = tp%dtheta
        tp%q0(2, :, :) = tp%dphi
        tp%q0(3:nq, :, :) = carried
        allocate (tp%ib1(nt, np), tp%w1(ns_tp + 1, nt, np))
        allocate (tp%p1(nq, ns_tp + 1, nt, np))
        do j = 1, np
            call moving_interp(tp%q0(:, :, j), 1, tv, ns_tp, tp%q1(:, :, j), &
                               tp%ib1(:, j), tp%w1(:, :, j), tp%p1(:, :, :, j))
        end do
        allocate (tp%ib2(np, nt), tp%w2(ns_tp + 1, np, nt))
        allocate (tp%p2(nq, ns_tp + 1, np, nt))
        do i = 1, nt
            call moving_interp(tp%q1(:, i, :), 2, pv, ns_tp, tp%q2(:, i, :), &
                               tp%ib2(:, i), tp%w2(:, :, i), tp%p2(:, :, :, i))
        end do
    end subroutine surface_forward

    !> Tangent of the surface transform. carried(nc, nt, np) are the carried
    !> quantities on the VMEC grid; outputs on the Boozer grid are
    !> q2(1:2+nc, nt, np) = (Delta theta, Delta varphi, carried).
    subroutine boozer_surface_jvp(lin, aiota, alam, bt, bp, carried, &
                                  aiota_d, alam_d, bt_d, bp_d, carried_d, &
                                  bcov_t, bcov_p, dtheta_V, dphi_V, q2, &
                                  bcov_t_d, bcov_p_d, dtheta_V_d, dphi_V_d, q2_d)
        type(boozer_surface_lin_t), intent(in) :: lin
        real(dp), intent(in) :: aiota, alam(:, :), bt(:, :), bp(:, :)
        real(dp), intent(in) :: carried(:, :, :)
        real(dp), intent(in) :: aiota_d, alam_d(:, :), bt_d(:, :), bp_d(:, :)
        real(dp), intent(in) :: carried_d(:, :, :)
        real(dp), intent(out) :: bcov_t, bcov_p, dtheta_V(:, :), dphi_V(:, :)
        real(dp), intent(out) :: q2(:, :, :)
        real(dp), intent(out) :: bcov_t_d, bcov_p_d, dtheta_V_d(:, :)
        real(dp), intent(out) :: dphi_V_d(:, :), q2_d(:, :, :)
        type(surface_tape_t) :: tp
        integer :: nt, np, nq, ns_tp, i, j
        real(dp) :: cells, dj_d, gbeg_d
        real(dp), allocatable :: raw_d(:, :), q0_d(:, :, :), q1_d(:, :, :)

        call surface_forward(lin, aiota, alam, bt, bp, carried, tp)
        nt = lin%grid%n_theta
        np = lin%grid%n_phi
        ns_tp = lin%grid%ns_tp
        nq = size(carried, 1) + 2
        cells = real((nt - 1)*(np - 1), dp)

        bcov_t = tp%bcov_t
        bcov_p = tp%bcov_p
        dtheta_V = tp%dtheta
        dphi_V = tp%dphi
        q2 = tp%q2

        bcov_t_d = sum(bt_d(2:nt, 2:np))/cells
        bcov_p_d = sum(bp_d(2:nt, 2:np))/cells
        dj_d = -tp%denomjac**2*(aiota_d*tp%bcov_t + aiota*bcov_t_d + bcov_p_d)
        gbeg_d = (bcov_t_d*tp%denomjac + tp%bcov_t*dj_d)*alam(1, 1) &
                 + tp%bcov_t*tp%denomjac*alam_d(1, 1)

        allocate (raw_d(nt, np))
        raw_d(:, 1) = matmul(lin%k_theta, bt_d(:, 1))
        do i = 1, nt
            raw_d(i, :) = raw_d(i, 1) + matmul(lin%k_phi, bp_d(i, :))
        end do
        dphi_V_d = dj_d*(tp%raw - tp%bcov_t*alam) &
                   + tp%denomjac*(raw_d - bcov_t_d*alam - tp%bcov_t*alam_d) + gbeg_d
        dtheta_V_d = aiota_d*tp%dphi + aiota*dphi_V_d + alam_d

        allocate (q0_d(nq, nt, np), q1_d(nq, nt, np))
        q0_d(1, :, :) = dtheta_V_d
        q0_d(2, :, :) = dphi_V_d
        q0_d(3:nq, :, :) = carried_d
        do j = 1, np
            call moving_interp_jvp(q0_d(:, :, j), 1, ns_tp, tp%ib1(:, j), &
                                   tp%w1(:, :, j), tp%p1(:, :, :, j), q1_d(:, :, j))
        end do
        do i = 1, nt
            call moving_interp_jvp(q1_d(:, i, :), 2, ns_tp, tp%ib2(:, i), &
                                   tp%w2(:, :, i), tp%p2(:, :, :, i), q2_d(:, i, :))
        end do
    end subroutine boozer_surface_jvp

    !> Adjoint of the surface transform: given cotangents of the outputs of
    !> boozer_surface_jvp, accumulate (+=) cotangents of the inputs.
    subroutine boozer_surface_vjp(lin, aiota, alam, bt, bp, carried, &
                                  bcov_t_b, bcov_p_b, dtheta_V_b, dphi_V_b, q2_b, &
                                  aiota_b, alam_b, bt_b, bp_b, carried_b)
        type(boozer_surface_lin_t), intent(in) :: lin
        real(dp), intent(in) :: aiota, alam(:, :), bt(:, :), bp(:, :)
        real(dp), intent(in) :: carried(:, :, :)
        real(dp), intent(in) :: bcov_t_b, bcov_p_b, dtheta_V_b(:, :)
        real(dp), intent(in) :: dphi_V_b(:, :), q2_b(:, :, :)
        real(dp), intent(inout) :: aiota_b, alam_b(:, :), bt_b(:, :), bp_b(:, :)
        real(dp), intent(inout) :: carried_b(:, :, :)
        type(surface_tape_t) :: tp
        integer :: nt, np, nq, ns_tp, i, j
        real(dp) :: cells, dj_b, gbeg_b, bt_mean_b, bp_mean_b, bt_lam_b
        real(dp), allocatable :: q0_b(:, :, :), q1_b(:, :, :), raw_b(:, :)
        real(dp), allocatable :: dphi_b(:, :), dtheta_b(:, :), raw1_b(:)

        call surface_forward(lin, aiota, alam, bt, bp, carried, tp)
        nt = lin%grid%n_theta
        np = lin%grid%n_phi
        ns_tp = lin%grid%ns_tp
        nq = size(carried, 1) + 2
        cells = real((nt - 1)*(np - 1), dp)

        allocate (q1_b(nq, nt, np), q0_b(nq, nt, np))
        q1_b = 0.0_dp
        q0_b = 0.0_dp
        do i = 1, nt
            call moving_interp_vjp(q2_b(:, i, :), 2, ns_tp, tp%ib2(:, i), &
                                   tp%w2(:, :, i), tp%p2(:, :, :, i), q1_b(:, i, :))
        end do
        do j = 1, np
            call moving_interp_vjp(q1_b(:, :, j), 1, ns_tp, tp%ib1(:, j), &
                                   tp%w1(:, :, j), tp%p1(:, :, :, j), q0_b(:, :, j))
        end do
        allocate (dtheta_b(nt, np), dphi_b(nt, np), raw_b(nt, np))
        carried_b = carried_b + q0_b(3:nq, :, :)
        dtheta_b = q0_b(1, :, :) + dtheta_V_b
        dphi_b = q0_b(2, :, :) + dphi_V_b

        aiota_b = aiota_b + sum(dtheta_b*tp%dphi)
        dphi_b = dphi_b + aiota*dtheta_b
        alam_b = alam_b + dtheta_b

        dj_b = sum(dphi_b*(tp%raw - tp%bcov_t*alam))
        gbeg_b = sum(dphi_b)
        raw_b = tp%denomjac*dphi_b
        alam_b = alam_b - tp%denomjac*tp%bcov_t*dphi_b
        bt_lam_b = -tp%denomjac*sum(dphi_b*alam)

        allocate (raw1_b(nt))
        do i = 1, nt
            raw1_b(i) = sum(raw_b(i, :))
            bp_b(i, :) = bp_b(i, :) + matmul(raw_b(i, :), lin%k_phi)
        end do
        bt_b(:, 1) = bt_b(:, 1) + matmul(raw1_b, lin%k_theta)

        bt_mean_b = bcov_t_b + bt_lam_b + gbeg_b*tp%denomjac*alam(1, 1)
        dj_b = dj_b + gbeg_b*tp%bcov_t*alam(1, 1)
        alam_b(1, 1) = alam_b(1, 1) + gbeg_b*tp%bcov_t*tp%denomjac

        aiota_b = aiota_b - dj_b*tp%denomjac**2*tp%bcov_t
        bt_mean_b = bt_mean_b - dj_b*tp%denomjac**2*aiota
        bp_mean_b = bcov_p_b - dj_b*tp%denomjac**2

        bt_b(2:nt, 2:np) = bt_b(2:nt, 2:np) + bt_mean_b/cells
        bp_b(2:nt, 2:np) = bp_b(2:nt, 2:np) + bp_mean_b/cells
    end subroutine boozer_surface_vjp

end module boozer_surface
