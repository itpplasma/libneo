module cylindrical_cartesian

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: cyl_to_cart, cart_to_cyl

contains

    pure subroutine cyl_to_cart(xfrom, xto, dxto_dxfrom)
        real(dp), intent(in) :: xfrom(3)
        real(dp), intent(out) :: xto(3)
        real(dp), intent(out), optional :: dxto_dxfrom(3,3)
        real(dp) :: R, phi
        real(dp) :: cos_phi, sin_phi

        R = xfrom(1)
        phi = xfrom(2)

        cos_phi = cos(phi)
        sin_phi = sin(phi)

        xto(1) = R * cos_phi
        xto(2) = R * sin_phi
        xto(3) = xfrom(3)

        if (present(dxto_dxfrom)) then
            dxto_dxfrom(1,1) = cos_phi
            dxto_dxfrom(1,2) = -R * sin_phi
            dxto_dxfrom(1,3) = 0.0_dp
            dxto_dxfrom(2,1) = sin_phi
            dxto_dxfrom(2,2) = R * cos_phi
            dxto_dxfrom(2,3) = 0.0_dp
            dxto_dxfrom(3,1) = 0.0_dp
            dxto_dxfrom(3,2) = 0.0_dp
            dxto_dxfrom(3,3) = 1.0_dp
        end if
    end subroutine cyl_to_cart

    pure subroutine cart_to_cyl(xfrom, xto, dxto_dxfrom)
        real(dp), intent(in) :: xfrom(3)
        real(dp), intent(out) :: xto(3)
        real(dp), intent(out), optional :: dxto_dxfrom(3,3)
        real(dp) :: x, y, z, r, phi

        x = xfrom(1)
        y = xfrom(2)
        z = xfrom(3)

        r = sqrt(x * x + y * y)
        phi = atan2(y, x)

        xto(1) = r
        xto(2) = phi
        xto(3) = z

        if (present(dxto_dxfrom)) then
            if (r > 0.0_dp) then
                dxto_dxfrom(1,1) = x / r
                dxto_dxfrom(1,2) = y / r
            else
                dxto_dxfrom(1,1) = 1.0_dp
                dxto_dxfrom(1,2) = 0.0_dp
            end if
            dxto_dxfrom(1,3) = 0.0_dp

            if (r > 0.0_dp) then
                dxto_dxfrom(2,1) = -y / (r * r)
                dxto_dxfrom(2,2) = x / (r * r)
            else
                dxto_dxfrom(2,1) = 0.0_dp
                dxto_dxfrom(2,2) = 0.0_dp
            end if
            dxto_dxfrom(2,3) = 0.0_dp

            dxto_dxfrom(3,1) = 0.0_dp
            dxto_dxfrom(3,2) = 0.0_dp
            dxto_dxfrom(3,3) = 1.0_dp
        end if
    end subroutine cart_to_cyl

end module cylindrical_cartesian

