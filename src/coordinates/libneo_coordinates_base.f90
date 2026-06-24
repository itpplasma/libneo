module libneo_coordinates_base
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_coordinate_conventions
    use math_constants, only: TWOPI
    use fortnum_status, only: fortnum_status_t
    use fortnum_multiroot, only: multiroot_hybrid
    implicit none
    private

    public :: coordinate_system_t
    public :: UNKNOWN, CYL, VMEC, BOOZER
    public :: RHO_TOR, RHO_POL, PSI_TOR_NORM, PSI_POL_NORM

    public :: chartmap_netcdf_spec, chartmap_from_cyl_ierr_spec
    public :: chartmap_from_cyl_ok, chartmap_from_cyl_err_max_iter
    public :: chartmap_from_cyl_err_singular, chartmap_from_cyl_err_out_of_bounds
    public :: chartmap_from_cyl_err_invalid
    public :: refcoords_file_unknown, refcoords_file_chartmap, refcoords_file_vmec_wout

    ! Geometric status of the Cartesian inverse cs%invert_cart, shared by the base
    ! generic inverse and every chart that overrides it (the enum lives on the base so
    ! a base-level inverse and the chartmap specialization report identical codes).
    ! These describe the geometry of the located root ONLY; they carry no
    ! confinement-loss meaning. A caller that needs a loss decision keys it on its own
    ! guiding-centre test, not on this status.
    !   LOCATED      converged interior root, residual at chart tolerance.
    !   CLAMPED_EDGE solver pinned rho to 1 and could not reduce the residual below
    !                tolerance, but the residual is still a small fraction of a radial
    !                cell: the target sits essentially on the last surface.
    !   OUTSIDE      best root has rho clamped to 1 with a residual that is a large
    !                fraction of a radial cell: the target is past the last surface.
    !   NO_ROOT      singular Jacobian or no convergence from any seed.
    integer, parameter, public :: CHARTMAP_LOCATED = 0
    integer, parameter, public :: CHARTMAP_CLAMPED_EDGE = 1
    integer, parameter, public :: CHARTMAP_OUTSIDE = 2
    integer, parameter, public :: CHARTMAP_NO_ROOT = 3

    ! A located root has residual below an absolute tolerance OR below this fraction of
    ! a radial cell |dx/drho| (scale-correct near the axis and on reactor-size charts
    ! where the absolute floor is meaningless). Shared by the generic inverse and the
    ! chartmap override so both classify a root the same way.
    real(dp), parameter, public :: coord_invert_edge_frac = 0.05_dp
    real(dp), parameter, public :: coord_invert_accept_tol = 1.0e-6_dp

    character(len=*), parameter :: chartmap_netcdf_spec = &
        "Chartmap NetCDF conventions:"//new_line('a')// &
        "- Dimensions: rho, theta, zeta"//new_line('a')// &
        "- Variables (required):"//new_line('a')// &
        "  rho(rho), theta(theta)"//new_line('a')// &
        "  zeta(zeta)"//new_line('a')// &
        "  x(zeta,theta,rho)"//new_line('a')// &
        "  y(zeta,theta,rho)"//new_line('a')// &
        "  z(zeta,theta,rho)"//new_line('a')// &
        "- Optional: num_field_periods (integer >= 1)"//new_line('a')// &
        "- Required: zeta_convention (global attribute)"//new_line('a')// &
        "  allowed: cyl, vmec, boozer"//new_line('a')// &
        "- Ranges:"//new_line('a')// &
        "  rho in [0,1]"//new_line('a')// &
        "  theta in [0,2pi)"//new_line('a')// &
        "  zeta in [0,2pi/num_field_periods)"//new_line('a')// &
        "periodic dims exclude endpoint"//new_line('a')// &
        "- Storage order:"//new_line('a')// &
        "  file dims (zeta,theta,rho)"//new_line('a')// &
        "Fortran reads x(rho,theta,zeta)"//new_line('a')// &
        "- Units: x,y,z in cm"//new_line('a')// &
        "- Contract:"//new_line('a')// &
        "  cyl/vmec require atan2(y,x)=zeta"//new_line('a')// &
        "  boozer stores true Cartesian geometry; zeta is a chart parameter"

    integer, parameter :: chartmap_from_cyl_ok = 0
    integer, parameter :: chartmap_from_cyl_err_max_iter = 1
    integer, parameter :: chartmap_from_cyl_err_singular = 2
    integer, parameter :: chartmap_from_cyl_err_out_of_bounds = 3
    integer, parameter :: chartmap_from_cyl_err_invalid = 4

    integer, parameter :: refcoords_file_unknown = 0
    integer, parameter :: refcoords_file_chartmap = 1
    integer, parameter :: refcoords_file_vmec_wout = 2

    character(len=*), parameter :: chartmap_from_cyl_ierr_spec = &
        "chartmap_from_cyl ierr codes:"//new_line('a')// &
        "0 ok"//new_line('a')// &
        "1 max iterations or no progress"//new_line('a')// &
        "2 singular normal equations"//new_line('a')// &
        "3 step out of bounds for rho"//new_line('a')// &
        "4 invalid mapping slice"

    type, abstract :: coordinate_system_t
    contains
        procedure(evaluate_cart_if), deferred :: evaluate_cart
        procedure(evaluate_cyl_if), deferred :: evaluate_cyl
        procedure(covariant_basis_if), deferred :: covariant_basis
        procedure(metric_tensor_if), deferred :: metric_tensor
        procedure(from_cyl_if), deferred :: from_cyl
        ! christoffel defaults to the finite-difference fallback so coordinate
        ! systems without an analytic form (e.g. flat Cartesian) need not override
        ! it; analytic charts (VMEC, geoflux, chartmap) override as before.
        procedure :: christoffel => coordinate_system_christoffel_default
        procedure :: cov_to_cart => coordinate_system_cov_to_cart
        procedure :: ctr_to_cart => coordinate_system_ctr_to_cart
        procedure :: christoffel_fd => coordinate_system_christoffel_fd
        ! Generic warm Cartesian inverse on the polymorphic evaluate_cart /
        ! covariant_basis (so any scaled override is used automatically). Charts that
        ! need an edge-bracketed specialization (chartmap) override this; VMEC and any
        ! other chart inherit the default below.
        procedure :: invert_cart => coordinate_system_invert_cart_default
    end type coordinate_system_t

    abstract interface
        subroutine evaluate_cart_if(self, u, x)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: x(3)
        end subroutine

        subroutine evaluate_cyl_if(self, u, x)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: x(3)
        end subroutine

        subroutine covariant_basis_if(self, u, e_cov)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: e_cov(3, 3)
        end subroutine

        subroutine metric_tensor_if(self, u, g, ginv, sqrtg)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: g(3, 3), ginv(3, 3), sqrtg
        end subroutine

        subroutine christoffel_if(self, u, Gamma)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: u(3)
            real(dp), intent(out) :: Gamma(3, 3, 3)
        end subroutine

        subroutine from_cyl_if(self, xcyl, u, ierr)
            import :: coordinate_system_t, dp
            class(coordinate_system_t), intent(in) :: self
            real(dp), intent(in) :: xcyl(3)
            real(dp), intent(out) :: u(3)
            integer, intent(out) :: ierr
        end subroutine
    end interface

    ! Read-only context threaded through the fortnum callback for the generic inverse:
    ! the polymorphic chart and the Cartesian target. The solver never inspects it; it
    ! passes it back to the callback unchanged. The pointer is polymorphic so a scaled
    ! subclass dispatches its own evaluate_cart / covariant_basis.
    type :: coord_invert_ctx_t
        class(coordinate_system_t), pointer :: self => null()
        real(dp) :: x_target(3) = 0.0_dp
    end type coord_invert_ctx_t

contains

    subroutine coordinate_system_cov_to_cart(self, u, v_cov, v_cart)
        class(coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(in) :: v_cov(3)
        real(dp), intent(out) :: v_cart(3)

        real(dp) :: e_cov(3, 3), g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: v_ctr(3)

        call self%metric_tensor(u, g, ginv, sqrtg)
        v_ctr = matmul(ginv, v_cov)
        call self%covariant_basis(u, e_cov)
        v_cart = matmul(e_cov, v_ctr)
    end subroutine coordinate_system_cov_to_cart

    subroutine coordinate_system_ctr_to_cart(self, u, v_ctr, v_cart)
        class(coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(in) :: v_ctr(3)
        real(dp), intent(out) :: v_cart(3)

        real(dp) :: e_cov(3, 3)

        call self%covariant_basis(u, e_cov)
        v_cart = matmul(e_cov, v_ctr)
    end subroutine coordinate_system_ctr_to_cart

    subroutine coordinate_system_christoffel_fd(self, u, Gamma, h)
        !> Christoffel symbols of the second kind from central differences of the
        !> coordinate system's own metric_tensor. Self-consistent fallback for
        !> charts without a closed-form metric derivative.
        !> Gamma(i,m,n) = Gamma^i_{mn} = 0.5 g^{il}(d_m g_{ln}+d_n g_{lm}-d_l g_{mn}).
        class(coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: Gamma(3, 3, 3)
        real(dp), intent(in), optional :: h

        real(dp) :: g(3, 3), ginv(3, 3), sqrtg
        real(dp) :: gp(3, 3), gm(3, 3), dummy_inv(3, 3), dummy_det
        real(dp) :: dg(3, 3, 3), up(3), um(3), hstep
        integer :: i, m, n, l, k

        hstep = 1.0e-5_dp
        if (present(h)) hstep = h

        call self%metric_tensor(u, g, ginv, sqrtg)

        do k = 1, 3
            up = u
            um = u
            up(k) = u(k) + hstep
            um(k) = u(k) - hstep
            call self%metric_tensor(up, gp, dummy_inv, dummy_det)
            call self%metric_tensor(um, gm, dummy_inv, dummy_det)
            dg(:, :, k) = (gp - gm)/(2.0_dp*hstep)
        end do

        Gamma = 0.0_dp
        do i = 1, 3
            do m = 1, 3
                do n = 1, 3
                    do l = 1, 3
                        Gamma(i, m, n) = Gamma(i, m, n) + 0.5_dp*ginv(i, l)* &
                            (dg(l, n, m) + dg(l, m, n) - dg(m, n, l))
                    end do
                end do
            end do
        end do
    end subroutine coordinate_system_christoffel_fd

    subroutine coordinate_system_christoffel_default(self, u, Gamma)
        !> Default christoffel binding: the finite-difference fallback. Matches the
        !> christoffel_if signature so analytic charts override it unchanged.
        class(coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp), intent(out) :: Gamma(3, 3, 3)

        call self%christoffel_fd(u, Gamma)
    end subroutine coordinate_system_christoffel_default

    ! Generic warm Cartesian inverse: from a carried guess u_guess locate the logical
    ! u = (rho, theta, phi) with self%evaluate_cart(u) = x. A damped Newton (fortnum
    ! multiroot_hybrid) runs in the pseudo-Cartesian chart w = (X, Y, phi) =
    ! (rho cos theta, rho sin theta, phi), which is regular through the magnetic axis:
    ! the polar column dx/dtheta ~ rho makes the polar Jacobian singular at rho = 0,
    ! while the (X, Y) Jacobian built by the chain rule stays finite. The callback
    ! evaluates the forward map and the analytic Jacobian through the POLYMORPHIC
    ! self%evaluate_cart / self%covariant_basis, so a scaled subclass inverts the scaled
    ! map without overriding this routine. geom_status is GEOMETRIC ONLY (located /
    ! clamped-edge / outside / no-root); it never encodes a confinement loss.
    subroutine coordinate_system_invert_cart_default(self, x, u_guess, u, geom_status)
        class(coordinate_system_t), intent(in), target :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(in) :: u_guess(3)
        real(dp), intent(out) :: u(3)
        integer, intent(out) :: geom_status

        type(coord_invert_ctx_t) :: ctx
        type(fortnum_status_t) :: status
        real(dp) :: w0(3), w(3), xc(3), res_norm, rscale

        ctx%self => self
        ctx%x_target = x

        w0(1) = u_guess(1)*cos(u_guess(2))
        w0(2) = u_guess(1)*sin(u_guess(2))
        w0(3) = u_guess(3)

        call multiroot_hybrid(coord_invert_fdf, 3, w0, w, status, &
                              xtol=1.0e-12_dp, ftol=1.0e-12_dp, max_iter=60, ctx=ctx)

        call coord_invert_w_to_u(w, u)
        u(1) = min(max(u(1), 0.0_dp), 1.0_dp)
        u(2) = modulo(u(2), TWOPI)

        call self%evaluate_cart(u, xc)
        res_norm = sqrt(sum((x - xc)**2))
        rscale = coordinate_system_radial_scale(self, u)
        geom_status = coordinate_system_classify_root(u(1), res_norm, rscale)
    end subroutine coordinate_system_invert_cart_default

    ! fortnum multiroot callback for the generic inverse: F(w) = x(w) - x_target with
    ! the analytic Jacobian dF/dw = pseudo-Cartesian Jacobian built from the polymorphic
    ! covariant_basis. A trial with rho = sqrt(X^2 + Y^2) > 1 is past the chart edge,
    ! where the forward map clamps rho and would report a spuriously small residual on
    ! the flat clamped region; penalize it (large outward F) so the line search
    ! backtracks instead of accepting the clamped-edge point.
    subroutine coord_invert_fdf(w, f, jac, ctx)
        real(dp), intent(in) :: w(:)
        real(dp), intent(out) :: f(:)
        real(dp), intent(out) :: jac(:, :)
        class(*), intent(in), optional :: ctx

        real(dp) :: u(3), e_cov(3, 3), jw(3, 3), xc(3), rho
        real(dp), parameter :: edge_penalty = 1.0e6_dp
        integer :: i

        f = 0.0_dp
        jac = 0.0_dp
        do i = 1, 3
            jac(i, i) = 1.0_dp
        end do
        if (.not. present(ctx)) return

        select type (ctx)
        type is (coord_invert_ctx_t)
            call coord_invert_w_to_u(w, u)
            rho = u(1)
            if (rho > 1.0_dp) then
                f = 0.0_dp
                f(1) = edge_penalty*(rho - 1.0_dp)
                return
            end if
            call ctx%self%evaluate_cart(u, xc)
            f = xc - ctx%x_target
            call ctx%self%covariant_basis(u, e_cov)
            call coord_invert_pseudocart_basis(u, e_cov, jw)
            jac = jw
        end select
    end subroutine coord_invert_fdf

    ! Pseudo-Cartesian near-axis chart basis jw(a,i) = dx_a/dw_i. w = (X, Y, phi) =
    ! (rho cos theta, rho sin theta, phi). dx/dX = cos theta dx/drho
    ! - (sin theta / rho) dx/dtheta, dx/dY = sin theta dx/drho
    ! + (cos theta / rho) dx/dtheta. The polar column dx/dtheta ~ rho cancels the
    ! 1/rho, so both columns stay finite at the axis.
    subroutine coord_invert_pseudocart_basis(u, e_cov, jw)
        real(dp), intent(in) :: u(3)
        real(dp), intent(in) :: e_cov(3, 3)
        real(dp), intent(out) :: jw(3, 3)
        real(dp) :: rho, cth, sth
        integer :: a

        rho = max(u(1), 1.0e-30_dp)
        cth = cos(u(2))
        sth = sin(u(2))
        do a = 1, 3
            jw(a, 1) = e_cov(a, 1)*cth - e_cov(a, 2)*(sth/rho)
            jw(a, 2) = e_cov(a, 1)*sth + e_cov(a, 2)*(cth/rho)
            jw(a, 3) = e_cov(a, 3)
        end do
    end subroutine coord_invert_pseudocart_basis

    ! Pseudo-Cartesian w = (X, Y, phi) -> logical u = (rho, theta, phi). rho >= 0 and
    ! the atan2 branch make the axis an ordinary point.
    pure subroutine coord_invert_w_to_u(w, u)
        real(dp), intent(in) :: w(3)
        real(dp), intent(out) :: u(3)
        u(1) = sqrt(w(1)**2 + w(2)**2)
        u(2) = atan2(w(2), w(1))
        u(3) = w(3)
    end subroutine coord_invert_w_to_u

    ! Length of one unit-rho radial step |dx/drho| = |e_cov(:,1)| through the
    ! polymorphic covariant_basis: the chart scale that sizes a residual.
    real(dp) function coordinate_system_radial_scale(self, u) result(s)
        class(coordinate_system_t), intent(in) :: self
        real(dp), intent(in) :: u(3)
        real(dp) :: e_cov(3, 3)
        call self%covariant_basis(u, e_cov)
        s = sqrt(e_cov(1, 1)**2 + e_cov(2, 1)**2 + e_cov(3, 1)**2)
    end function coordinate_system_radial_scale

    ! Classify a finished inverse from its radius and residual against the local radial
    ! cell. Interior at tolerance -> LOCATED; an edge-clamped root within a cell
    ! fraction -> CLAMPED_EDGE; a larger residual at rho = 1 -> OUTSIDE; otherwise
    ! NO_ROOT. The relative test keeps this scale-correct near the axis and on
    ! reactor-size charts.
    pure integer function coordinate_system_classify_root(rho, res_norm, rscale) &
        result(status)
        real(dp), intent(in) :: rho, res_norm, rscale
        logical :: located, at_edge

        located = res_norm < coord_invert_accept_tol .or. &
                  (rscale > 0.0_dp .and. res_norm < coord_invert_edge_frac*rscale)
        at_edge = rho >= 1.0_dp - coord_invert_edge_frac

        if (located .and. .not. at_edge) then
            status = CHARTMAP_LOCATED
        else if (located .and. at_edge) then
            status = CHARTMAP_CLAMPED_EDGE
        else if (at_edge) then
            status = CHARTMAP_OUTSIDE
        else
            status = CHARTMAP_NO_ROOT
        end if
    end function coordinate_system_classify_root

end module libneo_coordinates_base
