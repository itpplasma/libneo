module odeint_scale_floor_rhs
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: x_on = 1.0_dp
contains
    subroutine switch_on_rhs(x, y, dydx)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydx(:)

        dydx(1) = cos(x)
        dydx(2) = merge(1.0_dp, 0.0_dp, x > x_on)
    end subroutine switch_on_rhs
end module odeint_scale_floor_rhs

program test_odeint_scale_floor
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use odeint_allroutines_sub, only: odeint_allroutines, odeint_has_failed
    use odeint_scale_floor_rhs, only: x_on, switch_on_rhs
    implicit none

    real(dp), parameter :: eps = 1.0e-10_dp, x_end = 2.0_dp
    real(dp) :: y(2), y_ref(2)

    y_ref = [sin(x_end), x_end - x_on]

    y = [sin(x_on), 0.0_dp]
    call odeint_allroutines(y, 2, x_on, x_end, eps, switch_on_rhs)
    if (.not. odeint_has_failed()) then
        error stop 'zero-start component unexpectedly passed without scale floor'
    end if

    y = [sin(x_on), 0.0_dp]
    call odeint_allroutines(y, 2, x_on, x_end, eps, switch_on_rhs, &
                            scale_floor=[1.0_dp, 1.0_dp])
    if (odeint_has_failed()) error stop 'integration with scale floor failed'
    ! The jump at x_on costs one first-step error of about (b1/|e1|)*eps ~ 25*eps.
    if (any(abs(y - y_ref) > 100.0_dp*eps)) then
        write (*, *) 'error ', abs(y - y_ref)
        error stop 'integration with scale floor is inaccurate'
    end if
    write (*, *) 'odeint scale floor test passed'
end program test_odeint_scale_floor
