program test_polylag_5
use libneo_kinds, only : dp
use util_for_test, only: print_test, print_ok, print_fail
implicit none

call test_plag1d
call test_plag3d

contains

subroutine test_plag1d
use neo_polylag_5, only: find_node_index, plag1d
use libneo_util, only: linspace

real(dp), parameter :: tol = 1.0e-9_dp, x_min = 0.0_dp, x_max = 1.0_dp

integer :: nx
real(dp), dimension(:), allocatable :: x
real(dp) :: dx, x_test, f_test, df_test
integer, dimension(6) :: index
real(dp), dimension(6) :: xp, fp

call print_test("test_plag1d")

nx = int((x_max - x_min)/tol**(1.0_dp/5.0_dp))
allocate(x(nx))

call linspace(x_min, x_max, nx, x)
dx = (x(2) - x(1))

x_test = 0.5_dp
call find_node_index(x_test, x_min, dx, nx, index)

xp = x(index)
fp = f_1d(xp)

call plag1d(x_test, xp, fp, dx, f_test, df_test)

if (abs(f_test - f_1d(x_test)) > tol) then
    print *, "plag1d: ", f_test, " original: ", f_1d(x_test)
    call print_fail
    error stop
end if
if (abs(df_test - df_1d(x_test)) > tol) then
    print *, "plag1d derivative: ", df_test, " original: ", df_1d(x_test)
    call print_fail
    error stop
end if

call print_ok
end subroutine test_plag1d

elemental function f_1d(x) result(f)
    real(dp), intent(in) :: x
    real(dp) :: f

    f = sin(x)
end function f_1d

elemental function df_1d(x) result(df)
    real(dp), intent(in) :: x
    real(dp) :: df

    df = cos(x)
end function df_1d


subroutine test_plag3d
    use neo_polylag_5, only: find_node_index, get_plag3d
    use libneo_util, only: linspace

    real(dp), parameter :: tol = 1.0e-9_dp
    real(dp), parameter :: x_min = 0.0_dp, x_max = 1.0_dp
    real(dp), parameter :: y_min = 0.0_dp, y_max = 1.0_dp
    real(dp), parameter :: z_min = 0.0_dp, z_max = 1.0_dp

    integer :: nx, ny, nz
    real(dp), dimension(:), allocatable :: x, y, z
    real(dp) :: dx, dy, dz, x_test, y_test, z_test
    real(dp) :: f_test, dfdx_test, dfdy_test, dfdz_test
    integer, dimension(6) :: index
    real(dp), dimension(6) :: xp, yp, zp
    real(dp), dimension(6,6,6) :: fp
    integer :: i, j, k

    call print_test("test_plag3d")

    nx = int((x_max - x_min)/tol**(1.0_dp/5.0_dp))
    ny = int((y_max - y_min)/tol**(1.0_dp/5.0_dp))
    nz = int((z_max - z_min)/tol**(1.0_dp/5.0_dp))
    allocate(x(nx), y(ny), z(nz))

    call linspace(x_min, x_max, nx, x)
    call linspace(y_min, y_max, ny, y)
    call linspace(z_min, z_max, nz, z)

    dx = (x(2) - x(1))
    dy = (y(2) - y(1))
    dz = (z(2) - z(1))

    x_test = 0.5_dp
    y_test = 0.5_dp
    z_test = 0.5_dp
    call find_node_index(x_test, x_min, dx, nx, index)
    xp = x(index)
    call find_node_index(y_test, y_min, dy, ny, index)
    yp = y(index)
    call find_node_index(z_test, z_min, dz, nz, index)
    zp = z(index)

    do i = 1, 6
        do j = 1, 6
            do k = 1, 6
                fp(i,j,k) = f_3d(xp(i), yp(j), zp(k))
            end do
        end do
    end do

    call get_plag3d(x_test, y_test, z_test, xp, yp, zp, fp, dx, dy, dz, &
                f_test, dfdx_test, dfdy_test, dfdz_test)

    if (abs(f_test - f_3d(x_test, y_test, z_test)) > tol) then
        print *, "plag3d: ", f_test, " original: ", f_3d(x_test, y_test, z_test)
        call print_fail
        error stop
    end if
    if (abs(dfdx_test - df_3d(x_test, y_test, z_test, 'x')) > tol) then
        print *, "plag3d derivative x: ", dfdx_test, " original: ", df_3d(x_test, y_test, z_test, 'x')
        call print_fail
        error stop
    end if
    if (abs(dfdy_test - df_3d(x_test, y_test, z_test, 'y')) > tol) then
        print *, "plag3d derivative y: ", dfdy_test, " original: ", df_3d(x_test, y_test, z_test, 'y')
        call print_fail
        error stop
    end if
    if (abs(dfdz_test - df_3d(x_test, y_test, z_test, 'z')) > tol) then
        print *, "plag3d derivative z: ", dfdz_test, " original: ", df_3d(x_test, y_test, z_test, 'z')
        call print_fail
        error stop
    end if

    call print_ok
end subroutine test_plag3d

elemental function f_3d(x, y, z) result(f)
    real(dp), intent(in) :: x, y, z
    real(dp) :: f

    f = 0.5_dp * (sin(x) * cos(y) + exp(z-1))
end function f_3d

elemental function df_3d(x, y, z, variable) result(df)
    real(dp), intent(in) :: x, y ,z
    character(len=1), intent(in) :: variable
    real(dp) :: df

    select case (variable)
    case ('x')
        df = 0.5_dp * cos(x) * cos(y)
    case ('y')
        df = -0.5_dp * sin(x) * sin(y)
    case ('z')
        df = 0.5_dp * exp(z-1)
    case default
        df = 0.0_dp
    end select
end function df_3d

end program test_polylag_5
