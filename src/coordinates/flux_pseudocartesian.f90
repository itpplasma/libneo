module flux_pseudocartesian

    !> Pseudo-Cartesian near-axis chart for a flux radial coordinate s and a
    !> poloidal angle theta: X = sqrt(s) cos(theta), Y = sqrt(s) sin(theta), with
    !> the toroidal angle carried through unchanged. Flux coordinates are singular
    !> at the magnetic axis (s = 0): the poloidal angle is undefined and basis
    !> vectors diverge, so implicit orbit steps fail to converge through axis
    !> crossings. In (X, Y) the map back to (s, theta) is smooth (s = X^2 + Y^2,
    !> ds/dX = 2X, ...), so a solver that works in (X, Y) sees no singularity.
    !> See e.g. pseudo-Cartesian linearisation near the axis in guiding-centre
    !> orbit work (Pfefferle et al., arXiv:1412.5464).

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: flux_to_pseudocart, pseudocart_to_flux

contains

    pure subroutine flux_to_pseudocart(xfrom, xto, dxto_dxfrom)
        !> (s, theta, phi) -> (X, Y, phi). The radial column of the Jacobian
        !> carries the 1/(2 sqrt(s)) axis singularity of the forward map; it is
        !> set to zero at s = 0. Prefer pseudocart_to_flux for near-axis work.
        real(dp), intent(in) :: xfrom(3)
        real(dp), intent(out) :: xto(3)
        real(dp), intent(out), optional :: dxto_dxfrom(3, 3)
        real(dp) :: s, theta, rho, cth, sth

        s = xfrom(1)
        theta = xfrom(2)
        rho = sqrt(max(s, 0.0_dp))
        cth = cos(theta)
        sth = sin(theta)

        xto(1) = rho*cth
        xto(2) = rho*sth
        xto(3) = xfrom(3)

        if (present(dxto_dxfrom)) then
            dxto_dxfrom = 0.0_dp
            if (rho > 0.0_dp) then
                dxto_dxfrom(1, 1) = 0.5_dp*cth/rho
                dxto_dxfrom(2, 1) = 0.5_dp*sth/rho
            end if
            dxto_dxfrom(1, 2) = -rho*sth
            dxto_dxfrom(2, 2) = rho*cth
            dxto_dxfrom(3, 3) = 1.0_dp
        end if
    end subroutine flux_to_pseudocart

    pure subroutine pseudocart_to_flux(xfrom, xto, dxto_dxfrom)
        !> (X, Y, phi) -> (s, theta, phi). The map and its s-derivatives are
        !> smooth across the axis; only the theta-derivative is singular at the
        !> axis and is set to zero there.
        real(dp), intent(in) :: xfrom(3)
        real(dp), intent(out) :: xto(3)
        real(dp), intent(out), optional :: dxto_dxfrom(3, 3)
        real(dp) :: x, y, s

        x = xfrom(1)
        y = xfrom(2)
        s = x*x + y*y

        xto(1) = s
        xto(2) = atan2(y, x)
        xto(3) = xfrom(3)

        if (present(dxto_dxfrom)) then
            dxto_dxfrom = 0.0_dp
            dxto_dxfrom(1, 1) = 2.0_dp*x
            dxto_dxfrom(1, 2) = 2.0_dp*y
            if (s > 0.0_dp) then
                dxto_dxfrom(2, 1) = -y/s
                dxto_dxfrom(2, 2) = x/s
            end if
            dxto_dxfrom(3, 3) = 1.0_dp
        end if
    end subroutine pseudocart_to_flux

end module flux_pseudocartesian
