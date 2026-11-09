module neo_polylag_5
    use plag_coeff_sub, only: plag_coeff_uniform, lagrange_stencil
    use libneo_lagrange_point, only: lagrange_point_3d_n6, lagrange_value_3d_n6
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    integer, parameter :: mp = 6

contains

    subroutine find_node_index(u, u_min, du, nu, u_index)
! Stencil of mp nodes on the uniform mesh u_min + (i - 1)*du, i = 1..nu,
! starting at the cell that contains u.
        real(dp), intent(in) :: u, u_min, du
        integer, intent(in) :: nu
        integer, dimension(mp), intent(out) :: u_index

        call lagrange_stencil(int((u - u_min)/du) + 1, nu, u_index)
    end subroutine find_node_index

!> Lagrange interpolation on uniform (increasingly ordered) mesh 1D
    subroutine plag1d(x, xp, fp, dx, f, dfdx)
        real(dp), intent(in) :: x, dx
        real(dp), dimension(mp), intent(in) :: xp, fp
        real(dp), intent(out) :: f, dfdx

        real(dp) :: c(0:1, mp)

        call plag_coeff_uniform(mp, 1, x, xp(1), 1.0_dp/dx, c)
        f = dot_product(fp, c(0, :))
        dfdx = dot_product(fp, c(1, :))
    end subroutine plag1d

!> Lagrange interpolation on uniform (increasingly ordered) mesh 3D
    subroutine get_plag3d(x, y, z, xp, yp, zp, fp, dx, dy, dz, &
                          f, dfdx, dfdy, dfdz)
        real(dp), intent(in) :: x, y, z, dx, dy, dz
        real(dp), dimension(mp), intent(in) :: xp, yp, zp
        real(dp), dimension(mp, mp, mp), intent(in) :: fp
        real(dp), intent(out) :: f, dfdx, dfdy, dfdz

        real(dp) :: v(1), g(3, 1)

        call lagrange_point_3d_n6(fp, (x - xp(1))/dx, (y - yp(1))/dy, (z - zp(1))/dz, &
                                  1.0_dp/dx, 1.0_dp/dy, 1.0_dp/dz, v, g)
        f = v(1)
        dfdx = g(1, 1)
        dfdy = g(2, 1)
        dfdz = g(3, 1)
    end subroutine get_plag3d

    subroutine get_f_plag3d(x, y, z, xp, yp, zp, fp, dx, dy, dz, f)
        real(dp), intent(in) :: x, y, z, dx, dy, dz
        real(dp), dimension(mp), intent(in) :: xp, yp, zp
        real(dp), dimension(mp, mp, mp), intent(in) :: fp
        real(dp), intent(out) :: f

        real(dp) :: v(1)

        call lagrange_value_3d_n6(fp, (x - xp(1))/dx, (y - yp(1))/dy, (z - zp(1))/dz, &
                                  1.0_dp/dx, 1.0_dp/dy, 1.0_dp/dz, v)
        f = v(1)
    end subroutine get_f_plag3d

    subroutine get_df_plag3d(x, y, z, xp, yp, zp, fp, dx, dy, dz, dfdx, dfdy, dfdz)
        real(dp), intent(in) :: x, y, z, dx, dy, dz
        real(dp), dimension(mp), intent(in) :: xp, yp, zp
        real(dp), dimension(mp, mp, mp), intent(in) :: fp
        real(dp), intent(out) :: dfdx, dfdy, dfdz

        real(dp) :: f

        call get_plag3d(x, y, z, xp, yp, zp, fp, dx, dy, dz, f, dfdx, dfdy, dfdz)
    end subroutine get_df_plag3d

    function calc_coefs(u, up, du) result(coefs)
        real(dp), intent(in) :: u, du
        real(dp), dimension(mp), intent(in) :: up
        real(dp), dimension(mp) :: coefs

        real(dp) :: c(0:0, mp)

        call plag_coeff_uniform(mp, 0, u, up(1), 1.0_dp/du, c)
        coefs = c(0, :)
    end function calc_coefs

    function calc_derivative_coefs(u, up, du) result(coefs)
        real(dp), intent(in) :: u, du
        real(dp), dimension(mp), intent(in) :: up
        real(dp), dimension(mp) :: coefs

        real(dp) :: c(0:1, mp)

        call plag_coeff_uniform(mp, 1, u, up(1), 1.0_dp/du, c)
        coefs = c(1, :)
    end function calc_derivative_coefs

end module neo_polylag_5
