!> Cubic (mp = 4 point) Lagrange interpolation on uniform meshes. Legacy
!> interface kept for field_mesh3d, bdivfree and external callers; the weights
!> come from plag_coeff_uniform, the 3D evaluation from a generated kernel.
module polylag_3
    use libneo_kinds, only: dp
    use plag_coeff_sub, only: plag_coeff_uniform, lagrange_stencil
    use libneo_lagrange_point, only: lagrange_point_3d_n4
    implicit none

    integer, parameter :: mp = 4

contains

    !> Stencil around u on the uniform mesh umin + (i - 1)/dum1, i = 1..nup.
    subroutine indef(u, umin, dum1, nup, indu)
        real(dp), intent(in) :: u, umin, dum1
        integer, intent(in) :: nup
        integer, intent(out) :: indu(mp)

        call lagrange_stencil(int((u - umin)*dum1), nup, indu)
    end subroutine indef

    !> Stencil around mesh cell index.
    subroutine indsmp(index, nup, indu)
        integer, intent(in) :: index, nup
        integer, intent(out) :: indu(mp)

        call lagrange_stencil(index - 1, nup, indu)
    end subroutine indsmp

    subroutine plag1d(x, fp, dxm1, xp, polyl1d, p1x1d)
        real(dp), intent(in) :: x, fp(mp), dxm1, xp(mp)
        real(dp), intent(out) :: polyl1d, p1x1d

        real(dp) :: c(0:1, mp)

        call plag_coeff_uniform(mp, 1, x, xp(1), dxm1, c)
        polyl1d = sum(fp*c(0, :))
        p1x1d = sum(fp*c(1, :))
    end subroutine plag1d

    subroutine plag3d(x, y, z, fp, dxm1, dym1, dzm1, xp, yp, zp, &
                      polyl3d, poly1x, poly1y, poly1z)
        real(dp), intent(in) :: x, y, z, fp(mp, mp, mp), dxm1, dym1, dzm1
        real(dp), intent(in) :: xp(mp), yp(mp), zp(mp)
        real(dp), intent(out) :: polyl3d, poly1x, poly1y, poly1z

        real(dp) :: v(1), g(3, 1)

        call lagrange_point_3d_n4(fp, (x - xp(1))*dxm1, (y - yp(1))*dym1, &
                                  (z - zp(1))*dzm1, dxm1, dym1, dzm1, v, g)
        polyl3d = v(1)
        poly1x = g(1, 1)
        poly1y = g(2, 1)
        poly1z = g(3, 1)
    end subroutine plag3d

    subroutine coefs(u, up, dum1, cu)
        real(dp), intent(in) :: u, up(mp), dum1
        real(dp), intent(out) :: cu(mp)

        real(dp) :: c(0:0, mp)

        call plag_coeff_uniform(mp, 0, u, up(1), dum1, c)
        cu = c(0, :)
    end subroutine coefs

    subroutine coefs1(u, up, dum1, cu1)
        real(dp), intent(in) :: u, up(mp), dum1
        real(dp), intent(out) :: cu1(mp)

        real(dp) :: c(0:1, mp)

        call plag_coeff_uniform(mp, 1, u, up(1), dum1, c)
        cu1 = c(1, :)
    end subroutine coefs1

end module polylag_3
