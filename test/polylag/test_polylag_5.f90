program test_polylag_5
use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none

call test_plag1d

contains

subroutine test_plag1d
use neo_polylag_5, only: indef, plag1d
use util, only: linspace

integer, parameter :: nx = 10
real(dp), parameter :: tol = 1.0e-6_dp

real(dp), dimension(nx) :: x
real(dp) :: x_min, one_over_dx, x_test, f_test, df_test
integer, dimension(6) :: index
real(dp), dimension(6) :: xp, fp

call linspace(0.0_dp, 1.0_dp, nx, x)
x_min = minval(x)
one_over_dx = 1.0_dp / (x(2) - x(1))

x_test = 0.5_dp
call indef(x_test, x_min, one_over_dx, nx, index)

xp = x(index-2)
print *, "xp: ", xp
fp = example_function(xp)

call plag1d(x_test, fp, one_over_dx, xp, f_test, df_test)

if (abs(f_test - example_function(x_test)) > tol) then
    print *, "plag1d: ", f_test, " original: ", example_function(x_test)
    error stop
end if
if (abs(df_test - example_function_derivative(x_test)) > tol) then
    print *, "plag1d derivative: ", df_test, " original: ", example_function_derivative(x_test)
    error stop
end if
end subroutine test_plag1d

elemental function example_function(x) result(f)
    real(dp), intent(in) :: x
    real(dp) :: f

    f = sin(x)
end function example_function

elemental function example_function_derivative(x) result(df)
    real(dp), intent(in) :: x
    real(dp) :: df

    df = cos(x)
end function example_function_derivative
    
end program test_polylag_5