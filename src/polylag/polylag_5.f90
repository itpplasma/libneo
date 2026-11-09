module neo_polylag_5
    use plag_coeff_sub, only: plag_coeff_uniform, lagrange_stencil
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

        real(dp) :: cx(0:1, mp), cy(0:1, mp), cz(0:1, mp)
        real(dp) :: f0(mp, mp), f1(mp, mp), g00(mp), g10(mp), g01(mp)
        integer :: j, k

        call plag_coeff_uniform(mp, 1, x, xp(1), 1.0_dp/dx, cx)
        call plag_coeff_uniform(mp, 1, y, yp(1), 1.0_dp/dy, cy)
        call plag_coeff_uniform(mp, 1, z, zp(1), 1.0_dp/dz, cz)
        ! Contract one direction at a time: x, then y, then z.
        do k = 1, mp
            do j = 1, mp
                f0(j, k) = dot_product(fp(:, j, k), cx(0, :))
                f1(j, k) = dot_product(fp(:, j, k), cx(1, :))
            end do
        end do
        do k = 1, mp
            g00(k) = dot_product(f0(:, k), cy(0, :))
            g10(k) = dot_product(f1(:, k), cy(0, :))
            g01(k) = dot_product(f0(:, k), cy(1, :))
        end do
        f = dot_product(g00, cz(0, :))
        dfdx = dot_product(g10, cz(0, :))
        dfdy = dot_product(g01, cz(0, :))
        dfdz = dot_product(g00, cz(1, :))
    end subroutine get_plag3d

    subroutine get_f_plag3d(x, y, z, xp, yp, zp, fp, dx, dy, dz, f)
        real(dp), intent(in) :: x, y, z, dx, dy, dz
        real(dp), dimension(mp), intent(in) :: xp, yp, zp
        real(dp), dimension(mp, mp, mp), intent(in) :: fp
        real(dp), intent(out) :: f

        real(dp) :: cx(0:0, mp), cy(0:0, mp), cz(0:0, mp), f0(mp, mp), g(mp)
        integer :: j, k

        call plag_coeff_uniform(mp, 0, x, xp(1), 1.0_dp/dx, cx)
        call plag_coeff_uniform(mp, 0, y, yp(1), 1.0_dp/dy, cy)
        call plag_coeff_uniform(mp, 0, z, zp(1), 1.0_dp/dz, cz)
        do k = 1, mp
            do j = 1, mp
                f0(j, k) = dot_product(fp(:, j, k), cx(0, :))
            end do
            g(k) = dot_product(f0(:, k), cy(0, :))
        end do
        f = dot_product(g, cz(0, :))
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
