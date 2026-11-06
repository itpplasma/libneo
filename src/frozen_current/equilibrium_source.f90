!> Equilibrium sources of the frozen-current field.
!>
!> An equilibrium supplies the magnetic field B and the Jacobian at quadrature
!> points of its plasma volume V, given in flux coordinates (s, theta, zeta):
!> Gauss-Legendre in s on (0, 1), trapezoidal sums in theta on [0, 2 pi) and
!> in zeta over one field period [0, 2 pi/nfp), and on the boundary s = 1 with
!> its outward area vector.  Points carry cylindrical positions (R, phi, Z)
!> and physical cylindrical field components (B_R, B_phi, B_Z), in the
!> equilibrium's units (libneo's readers use cm and G).  The weights integrate
!> over one field period; the full torus is nfp periods.
module neo_equilibrium_source
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private

    public :: equilibrium_source_t, source_points_t, source_resolution_t
    public :: flux_quadrature

    real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp

    type :: source_resolution_t
        integer :: ns = 32 !! Gauss points in s
        integer :: ntheta = 64 !! trapezoidal points in theta
        integer :: nzeta = 32 !! trapezoidal points per field period
    end type source_resolution_t

    type :: source_points_t
        integer :: nfp = 1
        real(dp), allocatable :: x(:, :) !! (R, phi, Z)
        real(dp), allocatable :: b(:, :) !! (B_R, B_phi, B_Z)
        real(dp), allocatable :: sqg(:) !! oriented Jacobian of (s, theta, zeta)
        real(dp), allocatable :: w(:) !! volume: |sqg| ds dtheta dzeta
        real(dp), allocatable :: area(:, :) !! surface: outward n dA
    end type source_points_t

    type, abstract :: equilibrium_source_t
        integer :: nfp = 1
    contains
        procedure(sample_volume_i), deferred :: sample_volume
        procedure(sample_surface_i), deferred :: sample_surface
        procedure(boundary_curve_i), deferred :: boundary_curve
    end type equilibrium_source_t

    abstract interface
        subroutine sample_volume_i(self, res, pts)
            import :: equilibrium_source_t, source_resolution_t, source_points_t
            class(equilibrium_source_t), intent(in) :: self
            type(source_resolution_t), intent(in) :: res
            type(source_points_t), intent(out) :: pts
        end subroutine sample_volume_i

        subroutine sample_surface_i(self, res, pts)
            import :: equilibrium_source_t, source_resolution_t, source_points_t
            class(equilibrium_source_t), intent(in) :: self
            type(source_resolution_t), intent(in) :: res
            type(source_points_t), intent(out) :: pts
        end subroutine sample_surface_i

        !> Boundary s = 1 in the plane phi at the angles theta: (R, Z).
        subroutine boundary_curve_i(self, phi, theta, r, z)
            import :: equilibrium_source_t, dp
            class(equilibrium_source_t), intent(in) :: self
            real(dp), intent(in) :: phi, theta(:)
            real(dp), intent(out) :: r(:), z(:)
        end subroutine boundary_curve_i
    end interface

contains

    !> Nodes and weights of the flux-coordinate quadrature: s and ws from
    !> Gauss-Legendre on (0, 1), theta and zeta equispaced with the
    !> trapezoidal weights dtheta and dzeta of one field period.
    subroutine flux_quadrature(res, nfp, s, ws, theta, zeta, dtheta, dzeta)
        type(source_resolution_t), intent(in) :: res
        integer, intent(in) :: nfp
        real(dp), allocatable, intent(out) :: s(:), ws(:), theta(:), zeta(:)
        real(dp), intent(out) :: dtheta, dzeta
        integer :: i

        allocate (s(res%ns), ws(res%ns), theta(res%ntheta), zeta(res%nzeta))
        call gauss_legendre_ab(res%ns, 0.0_dp, 1.0_dp, s, ws)
        dtheta = twopi/res%ntheta
        dzeta = twopi/(real(nfp, dp)*res%nzeta)
        do i = 1, res%ntheta
            theta(i) = (i - 1)*dtheta
        end do
        do i = 1, res%nzeta
            zeta(i) = (i - 1)*dzeta
        end do
    end subroutine flux_quadrature

end module neo_equilibrium_source
