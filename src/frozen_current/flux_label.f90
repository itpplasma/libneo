!> Flux label continued beyond the plasma boundary.
!>
!> Inside the boundary s = 1 the label is the equilibrium's normalized
!> toroidal flux s, rho = sqrt(s).  Outside, in the plane phi = const,
!>     rho = 1 + d/a,    s = rho^2,
!> with d the distance to the boundary curve in that plane and
!> a = sqrt(A/pi) the effective minor radius of the cross-section of area A.
!> The label is continuous at the boundary, grows monotonically with the
!> distance, is defined at every point and needs no extrapolation of the
!> equilibrium's coordinates.  BEAMS3D (S_ARR), ASCOT5 (rho = sqrt(s), whose
!> hybrid field switches at rho = 1) and SIMPLE use this one convention.
module neo_flux_label
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_equilibrium_source, only: equilibrium_source_t
    implicit none
    private

    public :: boundary_section_t, boundary_section, continued_flux_label
    public :: inside_boundary

    real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp

    !> Boundary curve of one plane phi = const, closed polygon (R, Z).
    type :: boundary_section_t
        real(dp), allocatable :: r(:), z(:)
        real(dp) :: a = 0.0_dp !! effective minor radius sqrt(A/pi)
    end type boundary_section_t

contains

    subroutine boundary_section(src, phi, ntheta, sec)
        class(equilibrium_source_t), intent(in) :: src
        real(dp), intent(in) :: phi
        integer, intent(in) :: ntheta
        type(boundary_section_t), intent(out) :: sec
        real(dp) :: theta(ntheta), area
        integer :: j, k

        theta = [(twopi*(j - 1)/ntheta, j=1, ntheta)]
        allocate (sec%r(ntheta), sec%z(ntheta))
        call src%boundary_curve(phi, theta, sec%r, sec%z)
        area = 0.0_dp
        do j = 1, ntheta
            k = modulo(j, ntheta) + 1
            area = area + sec%r(j)*sec%z(k) - sec%r(k)*sec%z(j)
        end do
        sec%a = sqrt(0.5_dp*abs(area)/(0.5_dp*twopi))
    end subroutine boundary_section

    !> Crossing-number test of (r, z) against the closed boundary polygon.
    pure logical function inside_boundary(sec, r, z) result(inside)
        type(boundary_section_t), intent(in) :: sec
        real(dp), intent(in) :: r, z
        integer :: j, k, n

        n = size(sec%r)
        inside = .false.
        do j = 1, n
            k = modulo(j, n) + 1
            if ((sec%z(j) > z) .neqv. (sec%z(k) > z)) then
                if (r < sec%r(j) + (z - sec%z(j))*(sec%r(k) - sec%r(j))/ &
                    (sec%z(k) - sec%z(j))) inside = .not. inside
            end if
        end do
    end function inside_boundary

    !> Label s at (r, z) of the section's plane; s_inside is used inside the
    !> boundary, clipped to at most 1.
    pure real(dp) function continued_flux_label(sec, r, z, s_inside) result(s)
        type(boundary_section_t), intent(in) :: sec
        real(dp), intent(in) :: r, z, s_inside

        if (inside_boundary(sec, r, z)) then
            s = min(s_inside, 1.0_dp)
        else
            s = (1.0_dp + polygon_distance(sec, r, z)/sec%a)**2
        end if
    end function continued_flux_label

    pure real(dp) function polygon_distance(sec, r, z) result(d)
        type(boundary_section_t), intent(in) :: sec
        real(dp), intent(in) :: r, z
        real(dp) :: er, ez, t
        integer :: j, k, n

        n = size(sec%r)
        d = huge(1.0_dp)
        do j = 1, n
            k = modulo(j, n) + 1
            er = sec%r(k) - sec%r(j)
            ez = sec%z(k) - sec%z(j)
            t = ((r - sec%r(j))*er + (z - sec%z(j))*ez)/(er*er + ez*ez)
            t = max(0.0_dp, min(1.0_dp, t))
            d = min(d, hypot(r - sec%r(j) - t*er, z - sec%z(j) - t*ez))
        end do
    end function polygon_distance

end module neo_flux_label
