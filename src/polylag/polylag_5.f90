module neo_polylag_5
implicit none

integer, parameter :: dp = kind(1.0d0)

integer, parameter :: mp=6

contains

subroutine find_node_index(u,u_min,du,nu,u_index)
! defines interval for 1D interpolation on uniform mesh, normally
! looks for the central interval of stencil, but
! stops moving of stencil at the boundary (works for mp=4 only!)
! Input:
!    u - coordinate of a point (to be interpolated)
!    umin - minimal value of u
!    du = h - length of mesh interval
!    nu - total number of mesh points
! Output:
!    indu(mp) - relative index of stencil points
!
! the power 5 of polinomial is fixed strictly:
real(dp), intent(in) :: u, u_min, du
integer, intent(in) :: nu
integer, dimension(mp), intent(out) :: u_index

integer :: i

u_index(1) = int((u-u_min)/du)+1
if( u_index(1) .le. 0 ) u_index(1) = 1
u_index(mp) = u_index(1) + mp - 1
if( u_index(mp) .gt. nu ) then
    u_index(mp) = nu
    u_index(1) = u_index(mp) - mp + 1
endif
do i=2,mp-1
    u_index(i) = u_index(i-1) + 1
enddo
end subroutine find_node_index


!> Lagrange interpolation on uniform (increasingly ordered) mesh 1D
subroutine plag1d(x, xp, fp, dx, f, dfdx)
real(dp), intent(in) :: x, dx
real(dp), dimension(mp), intent(in) :: xp, fp
real(dp), intent(out) :: f, dfdx

real(dp), dimension(mp) :: x_coefs, dx_coefs
integer :: i

x_coefs = calc_coefs(x,xp,dx)
dx_coefs = calc_derivative_coefs(x,xp,dx)
f = 0.0_dp
dfdx = 0.0_dp
do i=1,mp
    f = f + fp(i) * x_coefs(i)
    dfdx = dfdx + fp(i) * dx_coefs(i)
enddo
end subroutine plag1d


!> Lagrange interpolation on uniform (increasingly ordered) mesh 3D
subroutine get_plag3d(x, y, z, xp, yp, zp, fp, dx, dy, dz, &
                  f, dfdx, dfdy, dfdz)
    real(dp), intent(in) :: x, y, z, dx, dy, dz
    real(dp), dimension(mp), intent(in) :: xp, yp, zp
    real(dp), dimension(mp,mp,mp), intent(in) :: fp
    real(dp), intent(out) :: f, dfdx, dfdy, dfdz

    call get_f_plag3d(x, y, z, xp, yp, zp, fp, dx, dy, dz, f)
    call get_df_plag3d(x, y, z, xp, yp, zp, fp, dx, dy, dz, dfdx, dfdy, dfdz)
end subroutine get_plag3d

subroutine get_f_plag3d(x, y, z, xp, yp, zp, fp, dx, dy, dz, f)
    real(dp), intent(in) :: x, y, z, dx, dy, dz
    real(dp), dimension(mp), intent(in) :: xp, yp, zp
    real(dp), dimension(mp,mp,mp), intent(in) :: fp
    real(dp), intent(out) :: f

    real(dp), dimension(mp) :: x_coefs, y_coefs, z_coefs
    integer :: i, j, k
    x_coefs = calc_coefs(x,xp,dx)
    y_coefs = calc_coefs(y,yp,dy)
    z_coefs = calc_coefs(z,zp,dz)
    f = 0.0_dp
    do k=1,mp
        do j=1,mp
            do i=1,mp
                f = f + fp(i,j,k) * x_coefs(i) * y_coefs(j) * z_coefs(k)
            enddo
        enddo
    enddo
end subroutine get_f_plag3d

subroutine get_df_plag3d(x, y, z, xp, yp, zp, fp, dx, dy, dz, dfdx, dfdy, dfdz)
    real(dp), intent(in) :: x, y, z, dx, dy, dz
    real(dp), dimension(mp), intent(in) :: xp, yp, zp
    real(dp), dimension(mp,mp,mp), intent(in) :: fp
    real(dp), intent(out) :: dfdx, dfdy, dfdz

    real(dp), dimension(mp) :: x_coefs, y_coefs, z_coefs
    real(dp), dimension(mp) :: dx_coefs, dy_coefs, dz_coefs
    integer :: i, j, k
    x_coefs = calc_coefs(x,xp,dx)
    y_coefs = calc_coefs(y,yp,dy)
    z_coefs = calc_coefs(z,zp,dz)
    dx_coefs = calc_derivative_coefs(x,xp,dx)
    dy_coefs = calc_derivative_coefs(y,yp,dy)
    dz_coefs = calc_derivative_coefs(z,zp,dz)
    dfdx = 0.0_dp
    dfdy = 0.0_dp
    dfdz = 0.0_dp
    do k=1,mp
        do j=1,mp
            do i=1,mp
                dfdx = dfdx + fp(i,j,k) * dx_coefs(i) * y_coefs(j) * z_coefs(k)
                dfdy = dfdy + fp(i,j,k) * x_coefs(i) * dy_coefs(j) * z_coefs(k)
                dfdz = dfdz + fp(i,j,k) * x_coefs(i) * y_coefs(j) * dz_coefs(k)
            enddo
        enddo
    enddo
end subroutine get_df_plag3d


function calc_coefs(u,up,du) result(coefs)
real(dp), intent(in) :: u, du
real(dp), dimension(mp), intent(in) :: up
real(dp), dimension(mp) :: coefs

real(dp) :: du5

du5 = du**5
coefs(1) = (u - up(2)) * (u - up(3)) * (u - up(4)) * &
           (u - up(5)) * (u - up(6)) / (-du5)/120.d0
coefs(2) = (u - up(1)) * (u - up(3)) * (u - up(4)) * &
           (u - up(5)) * (u - up(6)) / (du5)/24.d0
coefs(3) = (u - up(1)) * (u - up(2)) * (u - up(4)) * &
           (u - up(5)) * (u - up(6)) / (-du5)/12.d0
coefs(4) = (u - up(1)) * (u - up(2)) * (u - up(3)) * &
           (u - up(5)) * (u - up(6)) / (du5)/12.d0
coefs(5) = (u - up(1)) * (u - up(2)) * (u - up(3)) * &
           (u - up(4)) * (u - up(6)) / (-du5)/24.d0
coefs(6) = (u - up(1)) * (u - up(2)) * (u - up(3)) * &
           (u - up(4)) * (u - up(5)) / (du5)/120.d0
end function calc_coefs


function calc_derivative_coefs(u,up,du) result(coefs)
real(dp), intent(in) :: u, du
real(dp), dimension(mp), intent(in) :: up
real(dp), dimension(mp) :: coefs

real(dp) :: du5

du5 = du**5
coefs(1) = ((u - up(3)) * (u - up(4)) * &
            (u - up(5)) * (u - up(6)) + &
            (u - up(2)) * (u - up(4)) * &
            (u - up(5)) * (u - up(6)) + &
            (u - up(2)) * (u - up(3)) * &
            (u - up(5)) * (u - up(6)) + &
            (u - up(2)) * (u - up(3)) * &
            (u - up(4)) * (u - up(6)) + &
            (u - up(2)) * (u - up(3)) * &
            (u - up(4)) * (u - up(5)))/(-du5)/120.d0
coefs(2) = ((u - up(3)) * (u - up(4)) * &
            (u - up(5)) * (u - up(6)) + &
            (u - up(1)) * (u - up(4)) * &
            (u - up(5)) * (u - up(6)) + &
            (u - up(1)) * (u - up(3)) * &
            (u - up(5)) * (u - up(6)) + &
            (u - up(1)) * (u - up(3)) * &
            (u - up(4)) * (u - up(6)) + &
            (u - up(1)) * (u - up(3)) * &
            (u - up(4)) * (u - up(5)))/(du5)/24.d0
coefs(3) = ((u - up(2)) * (u - up(4)) * &
            (u - up(5)) * (u - up(6)) + &
            (u - up(1)) * (u - up(4)) * &
            (u - up(5)) * (u - up(6)) + &
            (u - up(1)) * (u - up(2)) * &
            (u - up(5)) * (u - up(6)) + &
            (u - up(1)) * (u - up(2)) * &
            (u - up(4)) * (u - up(6)) + &
            (u - up(1)) * (u - up(2)) * &
            (u - up(4)) * (u - up(5)))/(-du5)/12.d0
coefs(4) = ((u - up(2)) * (u - up(3)) * &
            (u - up(5)) * (u - up(6)) + &
            (u - up(1)) * (u - up(3)) * &
            (u - up(5)) * (u - up(6)) + &
            (u - up(1)) * (u - up(2)) * &
            (u - up(5)) * (u - up(6)) + &
            (u - up(1)) * (u - up(2)) * &
            (u - up(3)) * (u - up(6)) + &
            (u - up(1)) * (u - up(2)) * &
            (u - up(3)) * (u - up(5)))/(du5)/12.d0
coefs(5) = ((u - up(2)) * (u - up(3)) * &
            (u - up(4)) * (u - up(6)) + &
            (u - up(1)) * (u - up(3)) * &
            (u - up(4)) * (u - up(6)) + &
            (u - up(1)) * (u - up(2)) * &
            (u - up(4)) * (u - up(6)) + &
            (u - up(1)) * (u - up(2)) * &
            (u - up(3)) * (u - up(6)) + &
            (u - up(1)) * (u - up(2)) * &
            (u - up(3)) * (u - up(4)))/(-du5)/24.d0
coefs(6) = ((u - up(2)) * (u - up(3)) * &
            (u - up(4)) * (u - up(5)) + &
            (u - up(1)) * (u - up(3)) * &
            (u - up(4)) * (u - up(5)) + &
            (u - up(1)) * (u - up(2)) * &
            (u - up(4)) * (u - up(5)) + &
            (u - up(1)) * (u - up(2)) * &
            (u - up(3)) * (u - up(5)) + &
            (u - up(1)) * (u - up(2)) * &
            (u - up(3)) * (u - up(4)))/(du5)/120.d0
end function calc_derivative_coefs

end module neo_polylag_5
