!> Frozen-current magnetic field A = A_coils + A_p, B = B_coils + curl A_p.
!>
!> A_p is the vector potential of the plasma current held fixed, from the
!> weak problem of itpplasma/libneo#493: find A_p and the Coulomb-gauge
!> multiplier p with
!>     (curl A_p, curl v) + (grad p, v) = (B_p, curl v)_V + <n x B_p, v>_dV,
!>     (A_p, grad q) = 0,
!> where B_p = B - B_coils on the plasma volume V.  The right-hand side
!> equals (curl B_p, v)_V, so only the plasma current enters.  Each toroidal
!> mode n = k nfp, k = 0..kmax, is solved on a rectangle of the (R, Z) plane
!> with FortFEM's degree-p B-spline complex and its exact exterior
!> (fortfem_meridian_bspline_n), so the field holds from the axis through
!> the vacuum region inside the rectangle.  A source point of weight W at
!> toroidal angle phi enters mode n as W B exp(-i n phi) nfp/(2 pi), the
!> factor nfp summing the field periods.  The real field is
!>     A_p = Re a_0 + 2 Re sum_{k>0} a_{k nfp} exp(i k nfp phi).
!>
!> Units follow the inputs: cm and G with libneo's VMEC reader and coils.
!> The field_t interface takes Cartesian x and returns Cartesian A and B;
!> without coils it returns the plasma part alone.
module neo_frozen_current_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_meridian_bspline_n, only: meridian_bspline_t, &
                                meridian_bspline_build, meridian_bspline_exterior_build, &
                                     meridian_bspline_solve, meridian_bspline_load_curl, &
                                     meridian_bspline_load_value, meridian_bspline_eval, &
                                          meridian_bspline_field_energy, &
                                          meridian_bspline_eval_gradient
    use fortfem_meridian_magnetic_coupling, only: meridian_exterior_t
    use fortfem_meridian_feec_n, only: meridian_op_complex
    use fortsparse, only: csc_z_t
    use neo_field_base, only: field_t
    use neo_biotsavart, only: coils_t, compute_vector_potential, &
                              compute_magnetic_field, &
                              coil_field_jacobian => compute_field_jacobian
    use neo_equilibrium_source, only: source_points_t
    use neo_frozen_current_leaves, only: area_cross_field, cyl_to_cart
    implicit none
    private

    public :: frozen_current_patch_t, frozen_current_plasma_t
    public :: frozen_current_field_t, solve_frozen_current
    public :: frozen_current_plasma_cyl, frozen_current_plasma_cyl_gradient

    real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp

    !> Meridional rectangle, spline resolution and mode truncation.
    type :: frozen_current_patch_t
        real(dp) :: r0 = 0.0_dp, r1 = 0.0_dp, z0 = 0.0_dp, z1 = 0.0_dp
        integer :: nr = 64, nz = 64, degree = 3
        integer :: kmax = 8 !! modes n = k nfp, k = 0..kmax
        integer :: panels_per_span = 1
    end type frozen_current_patch_t

    !> Solved plasma potential: V1 coefficients per mode and diagnostics.
    type :: frozen_current_plasma_t
        type(meridian_bspline_t) :: sp
        integer :: nfp = 1, kmax = 0
        complex(dp), allocatable :: u(:, :) !! (n1, 0:kmax)
        real(dp), allocatable :: energy(:) !! (curl a_n, curl a_n)
        real(dp), allocatable :: gauge(:) !! |M1 grad p| / |f| per mode
    end type frozen_current_plasma_t

    type, extends(field_t) :: frozen_current_field_t
        type(coils_t) :: coils
        type(frozen_current_plasma_t) :: plasma
    contains
        procedure :: compute_afield
        procedure :: compute_bfield
        procedure :: compute_abfield
        procedure :: compute_field_jacobian
    end type frozen_current_field_t

contains

    !> Solve every mode for the plasma field B_p given at the volume points
    !> vol (weights w) and the boundary points surf (area vectors).
    subroutine solve_frozen_current(patch, vol, surf, plasma, status)
        type(frozen_current_patch_t), intent(in) :: patch
        type(source_points_t), intent(in) :: vol, surf
        type(frozen_current_plasma_t), intent(out) :: plasma
        integer, intent(out) :: status

        type(meridian_exterior_t) :: ex
        complex(dp), allocatable :: f(:), u(:), p(:), bv(:, :), gs(:, :)
        real(dp), allocatable :: xv(:, :), xs(:, :), ws(:), g(:, :)
        integer :: k, n, m

        call meridian_bspline_build(patch%r0, patch%r1, patch%z0, patch%z1, &
                                    patch%nr, patch%nz, patch%degree, plasma%sp, &
                                    status)
        if (status /= 0) return
        plasma%nfp = vol%nfp
        plasma%kmax = patch%kmax
        allocate (plasma%u(plasma%sp%n1, 0:patch%kmax))
        allocate (plasma%energy(0:patch%kmax), plasma%gauge(0:patch%kmax))

        allocate (xv(2, size(vol%x, 2)), xs(2, size(surf%x, 2)))
        xv = vol%x([1, 3], :)
        xs = surf%x([1, 3], :)
        allocate (ws(size(surf%x, 2)), g(3, size(surf%x, 2)))
        allocate (bv(3, size(vol%x, 2)), gs(3, size(surf%x, 2)), f(plasma%sp%n1))
        ws = 1.0_dp
        do m = 1, size(surf%x, 2)
            call area_cross_field(surf%area(1, m), surf%area(2, m), surf%area(3, m), &
                                  surf%b(1, m), surf%b(2, m), surf%b(3, m), g(1, m), &
                                  g(2, m), g(3, m))
        end do

        do k = 0, patch%kmax
            n = k*vol%nfp
            call mode_weights(vol%b, vol%x(2, :), n, vol%nfp, bv)
            call mode_weights(g, surf%x(2, :), n, vol%nfp, gs)
            f = (0.0_dp, 0.0_dp)
            call meridian_bspline_load_curl(plasma%sp, n, xv, vol%w, bv, f, status)
            if (status /= 0) return
            call meridian_bspline_load_value(plasma%sp, n, xs, ws, gs, f, status)
            if (status /= 0) return
            call meridian_bspline_exterior_build(plasma%sp, n, ex, status, &
                                                 patch%panels_per_span)
            if (status /= 0) return
            call meridian_bspline_solve(plasma%sp, ex, f, u, p, status)
            if (status /= 0) return
            plasma%u(:, k) = u
            call meridian_bspline_field_energy(plasma%sp, n, u, plasma%energy(k), &
                                               status)
            if (status /= 0) return
            plasma%gauge(k) = gauge_fraction(plasma%sp, n, p, f)
        end do
    end subroutine solve_frozen_current

    subroutine mode_weights(b, phi, n, nfp, bn)
        real(dp), intent(in) :: b(:, :), phi(:)
        integer, intent(in) :: n, nfp
        complex(dp), intent(out) :: bn(:, :)
        complex(dp) :: e
        integer :: m

        !$omp parallel do default(shared) private(m, e)
        do m = 1, size(b, 2)
            e = cmplx(cos(n*phi(m)), -sin(n*phi(m)), dp)*(nfp/twopi)
            bn(:, m) = b(:, m)*e
        end do
        !$omp end parallel do
    end subroutine mode_weights

    !> Size of the gauge term M1 grad_n p relative to the load f.  It vanishes
    !> for a divergence-free current tangential to the boundary.
    real(dp) function gauge_fraction(sp, n, p, f) result(r)
        type(meridian_bspline_t), intent(in) :: sp
        integer, intent(in) :: n
        complex(dp), intent(in) :: p(:), f(:)
        type(csc_z_t) :: grad
        complex(dp), allocatable :: gp(:), mgp(:)
        integer :: j, e, st

        r = 0.0_dp
        if (norm2(abs(f)) == 0.0_dp) return
        call meridian_op_complex(sp%grad, n, grad, st)
        if (st /= 0) error stop 'gauge_fraction: gradient operator'
        allocate (gp(grad%nrow), mgp(sp%mass(1)%nrow))
        gp = (0.0_dp, 0.0_dp)
        do j = 1, grad%ncol
            do e = grad%col_ptr(j), grad%col_ptr(j + 1) - 1
                gp(grad%row_idx(e)) = gp(grad%row_idx(e)) + grad%val(e)*p(j)
            end do
        end do
        mgp = (0.0_dp, 0.0_dp)
        do j = 1, sp%mass(1)%ncol
            do e = sp%mass(1)%col_ptr(j), sp%mass(1)%col_ptr(j + 1) - 1
                mgp(sp%mass(1)%row_idx(e)) = mgp(sp%mass(1)%row_idx(e)) + &
                                             sp%mass(1)%val(e)*gp(j)
            end do
        end do
        r = norm2(abs(mgp))/norm2(abs(f))
    end function gauge_fraction

    !> Plasma potential and field, physical cylindrical components, at the
    !> cylindrical point x = (R, phi, Z) inside the patch.
    subroutine frozen_current_plasma_cyl(plasma, x, a, b, status)
        type(frozen_current_plasma_t), intent(in) :: plasma
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: a(3), b(3)
        integer, intent(out) :: status

        complex(dp) :: an(3, 1), bn(3, 1), e
        real(dp) :: fac
        integer :: k, n

        a = 0.0_dp
        b = 0.0_dp
        do k = 0, plasma%kmax
            n = k*plasma%nfp
            call meridian_bspline_eval(plasma%sp, n, plasma%u(:, k), &
                                       reshape(x([1, 3]), [2, 1]), an, bn, status)
            if (status /= 0) return
            fac = merge(1.0_dp, 2.0_dp, k == 0)
            e = cmplx(cos(n*x(2)), sin(n*x(2)), dp)
            a = a + fac*real(an(:, 1)*e, dp)
            b = b + fac*real(bn(:, 1)*e, dp)
        end do
    end subroutine frozen_current_plasma_cyl

    !> Physical cylindrical components b = (B_R, B_phi, B_Z) of the plasma
    !> field and db(i, j), the derivative of b(i) along (R, phi, Z)(j).
    subroutine frozen_current_plasma_cyl_gradient(plasma, x, b, db, status)
        type(frozen_current_plasma_t), intent(in) :: plasma
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3), db(3, 3)
        integer, intent(out) :: status

        complex(dp) :: bn(3, 1), dbn(3, 2, 1), e
        real(dp) :: fac
        integer :: k, n

        b = 0.0_dp
        db = 0.0_dp
        do k = 0, plasma%kmax
            n = k*plasma%nfp
            call meridian_bspline_eval_gradient(plasma%sp, n, plasma%u(:, k), &
                                                reshape(x([1, 3]), [2, 1]), bn, dbn, &
                                                status)
            if (status /= 0) return
            fac = merge(1.0_dp, 2.0_dp, k == 0)
            e = cmplx(cos(n*x(2)), sin(n*x(2)), dp)
            b = b + fac*real(bn(:, 1)*e, dp)
            db(:, 1) = db(:, 1) + fac*real(dbn(:, 1, 1)*e, dp)
            db(:, 2) = db(:, 2) + fac*real(cmplx(0.0_dp, real(n, dp), dp)*bn(:, 1)*e, dp)
            db(:, 3) = db(:, 3) + fac*real(dbn(:, 2, 1)*e, dp)
        end do
    end subroutine frozen_current_plasma_cyl_gradient

    !> Cartesian B and dB(i, j) = dB_i/dx_j of plasma part and coils.
    subroutine compute_field_jacobian(self, x, B, dB)
        class(frozen_current_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: B(3), dB(3, 3)

        real(dp) :: b_cyl(3), db_cyl(3, 3), A(3), dA(3, 3), bc(3), dbc(3, 3)
        integer :: status

        call frozen_current_plasma_cyl_gradient(self%plasma, &
                                                [hypot(x(1), x(2)), atan2(x(2), x(1)), &
                                                 x(3)], b_cyl, db_cyl, status)
        if (status /= 0) error stop 'frozen_current_field: point outside the patch'
        call cylindrical_vector_jacobian(x(1), x(2), b_cyl(1), b_cyl(2), b_cyl(3), &
                                         db_cyl(1, 1), db_cyl(2, 1), db_cyl(3, 1), &
                                         db_cyl(1, 2), db_cyl(2, 2), db_cyl(3, 2), &
                                         db_cyl(1, 3), db_cyl(2, 3), db_cyl(3, 3), &
                                         B(1), B(2), B(3), dB(1, 1), dB(1, 2), &
                                         dB(1, 3), dB(2, 1), dB(2, 2), dB(2, 3), &
                                         dB(3, 1), dB(3, 2), dB(3, 3))
        if (.not. allocated(self%coils%x)) return
        call coil_field_jacobian(self%coils, x, A, dA, bc, dbc)
        B = B + bc
        dB = dB + dbc
    end subroutine compute_field_jacobian

    include '../generated/cylindrical_vector_jacobian_leaves.inc'

    subroutine compute_abfield(self, x, A, B)
        class(frozen_current_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: A(3), B(3)

        real(dp) :: r, c, s, ap(3), bp(3)
        integer :: status

        r = hypot(x(1), x(2))
        c = x(1)/r
        s = x(2)/r
        call frozen_current_plasma_cyl(self%plasma, [r, atan2(x(2), x(1)), x(3)], &
                                       ap, bp, status)
        if (status /= 0) error stop 'frozen_current_field: point outside the patch'
        call cyl_to_cart(c, s, ap(1), ap(2), A(1), A(2))
        call cyl_to_cart(c, s, bp(1), bp(2), B(1), B(2))
        A(3) = ap(3)
        B(3) = bp(3)
        if (.not. allocated(self%coils%x)) return
        A = A + compute_vector_potential(self%coils, x)
        B = B + compute_magnetic_field(self%coils, x)
    end subroutine compute_abfield

    subroutine compute_afield(self, x, A)
        class(frozen_current_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: A(3)
        real(dp) :: B(3)

        call self%compute_abfield(x, A, B)
    end subroutine compute_afield

    subroutine compute_bfield(self, x, B)
        class(frozen_current_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: B(3)
        real(dp) :: A(3)

        call self%compute_abfield(x, A, B)
    end subroutine compute_bfield

end module neo_frozen_current_field
