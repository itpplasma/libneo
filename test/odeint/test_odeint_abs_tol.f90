module odeint_abs_tol_rhs
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
contains
    ! y1' = cos(x); y2' = 1 for x > 0, else 0. At x = 0 the second component
    ! is zero with zero derivative, but its derivative is 1 at every later
    ! stage, so its embedded error estimate is |e1|*h > 0 for any h > 0.
    subroutine rhs_switch_on(x, y, dydx)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydx(:)
        dydx(1) = cos(x)
        dydx(2) = merge(1.0_dp, 0.0_dp, x > 0.0_dp)
    end subroutine rhs_switch_on

    ! Nested running integrals from zero, like field-line quadratures:
    ! y1 = sin(x), y2 = int y1 = 1 - cos(x), y3 = int y2 = x - sin(x).
    subroutine rhs_nested(x, y, dydx)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydx(:)
        dydx(1) = cos(x)
        dydx(2) = y(1)
        dydx(3) = y(2)
    end subroutine rhs_nested

    ! y' = y**2, y(0) = 1: y = 1/(1-x) blows up at x = 1.
    subroutine rhs_blowup(x, y, dydx)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydx(:)
        dydx(1) = y(1)**2
    end subroutine rhs_blowup

    ! y' = 1 for x <= 0.5, NaN beyond.
    subroutine rhs_nan_after_half(x, y, dydx)
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydx(:)
        dydx(1) = 1.0_dp
        if (x > 0.5_dp) dydx(1) = ieee_value(1.0_dp, ieee_quiet_nan)
    end subroutine rhs_nan_after_half

    subroutine rhs_constant(x, y, dydx)
        real(dp), intent(in) :: x, y(:)
        real(dp), intent(out) :: dydx(:)
        dydx = 1.0_dp
    end subroutine rhs_constant

    subroutine rhs_large_constant(x, y, dydx)
        real(dp), intent(in) :: x, y(:)
        real(dp), intent(out) :: dydx(:)
        dydx = 0.75_dp*huge(1.0_dp)
    end subroutine rhs_large_constant

    subroutine rhs_switch_on_ctx(x, y, dydx, context)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydx(:)
        class(*), intent(in) :: context
        call rhs_switch_on(x, y, dydx)
    end subroutine rhs_switch_on_ctx

    subroutine rhs_nan_ctx(x, y, dydx, context)
        real(dp), intent(in) :: x, y(:)
        real(dp), intent(out) :: dydx(:)
        class(*), intent(in) :: context
        call rhs_nan_after_half(x, y, dydx)
    end subroutine rhs_nan_ctx
end module odeint_abs_tol_rhs

program test_odeint_abs_tol
    use odeint_abs_tol_rhs
    use, intrinsic :: ieee_arithmetic, only: &
        ieee_value, ieee_quiet_nan, ieee_positive_inf
    use odeint_allroutines_sub, only: odeint_allroutines, &
        odeint_allroutines_checked, odeint_has_failed
    implicit none

    real(dp), parameter :: rtol = 1.0e-10_dp
    real(dp) :: y2(2), y3(3), y1(1), x, h, err
    integer :: i, ierr
    logical :: failed
    character(len=32) :: mode

    call get_command_argument(1, mode)
    if (trim(mode) == 'unchecked') then
        y2 = 0.0_dp
        call odeint_allroutines(y2, 2, 0.0_dp, 1.0_dp, rtol, rhs_switch_on)
        print *, 'FAIL: unchecked integration failure returned'
        stop
    end if
    if (trim(mode) == 'unchecked_context') then
        y2 = 0.0_dp
        call odeint_allroutines(y2, 2, 0, 0.0_dp, 1.0_dp, rtol, rhs_switch_on_ctx)
        print *, 'FAIL: unchecked context integration failure returned'
        stop
    end if

    if (trim(mode) == 'invalid_span') then
        y1 = 0.0_dp
        call odeint_allroutines_checked(y1, 1, ieee_value(1.0_dp, ieee_quiet_nan), &
            1.0_dp, rtol, rhs_constant, ierr)
        if (ierr /= 0) then
            print *, 'PASS: nonfinite bounds report failure'
        else
            print *, 'FAIL: nonfinite bounds accepted'
        end if
        stop
    end if
    if (trim(mode) == 'invalid_step') then
        y1 = 0.0_dp
        call odeint_allroutines_checked(y1, 1, 0.0_dp, 1.0_dp, rtol, &
            rhs_constant, ierr, &
            initial_stepsize=ieee_value(1.0_dp, ieee_quiet_nan))
        if (ierr /= 0) then
            print *, 'PASS: nonfinite step reports failure'
        else
            print *, 'FAIL: nonfinite step accepted'
        end if
        stop
    end if

    failed = .false.

    ! 1. Purely relative control cannot integrate the switch-on component.
    y2 = 0.0_dp
    call odeint_allroutines_checked(y2, 2, 0.0_dp, 1.0_dp, rtol, rhs_switch_on, ierr)
    call check(ierr == 1 .and. odeint_has_failed(), &
        'relative control reports failure on zero start component')

    ! 2. Mixed control integrates it to the analytic solution.
    y2 = 0.0_dp
    call odeint_allroutines_checked(y2, 2, 0.0_dp, 1.0_dp, rtol, rhs_switch_on, ierr, &
        atol=[rtol, rtol])
    err = max(abs(y2(1) - sin(1.0_dp)), abs(y2(2) - 1.0_dp))
    call check(ierr == 0 .and. .not. odeint_has_failed() .and. err < 1.0e-8_dp, &
        'mixed control integrates switch-on component')

    ! 3. Nested running integrals over many output intervals.
    y3 = 0.0_dp
    x = 0.0_dp
    h = 0.1_dp
    do i = 1, 100
        call odeint_allroutines_checked(y3, 3, x, x + h, rtol, rhs_nested, ierr, &
            atol=[rtol, rtol, rtol])
        if (ierr /= 0) exit
        x = x + h
    end do
    err = maxval(abs(y3 - [sin(x), 1.0_dp - cos(x), x - sin(x)]))
    call check(ierr == 0 .and. .not. odeint_has_failed() .and. err < 1.0e-8_dp, &
        'mixed control matches nested integrals')

    ! 4. Excluded component (atol = huge) does not break the others.
    y3 = 0.0_dp
    call odeint_allroutines_checked(y3, 3, 0.0_dp, 1.0_dp, rtol, rhs_nested, ierr, &
        atol=[rtol, rtol, huge(1.0_dp)])
    err = maxval(abs(y3(1:2) - [sin(1.0_dp), 1.0_dp - cos(1.0_dp)]))
    call check(ierr == 0 .and. .not. odeint_has_failed() .and. err < 1.0e-8_dp, &
        'excluded component')

    ! 5. Failure status with atol: a zero tolerance for the switch-on
    !    component is the purely relative test again.
    y2 = 0.0_dp
    call odeint_allroutines_checked(y2, 2, 0.0_dp, 1.0_dp, rtol, rhs_switch_on, ierr, &
        atol=[rtol, 0.0_dp])
    call check(ierr == 1 .and. odeint_has_failed(), 'zero atol reports failure')

    ! 5b. A non-finite right-hand side beyond x = 0.5 must not be accepted.
    y1 = 0.0_dp
    call odeint_allroutines_checked(y1, 1, 0.0_dp, 1.0_dp, rtol, &
        rhs_nan_after_half, ierr, &
        atol=[rtol])
    call check(odeint_has_failed() .and. abs(y1(1) - 0.5_dp) < 1.0e-8_dp, &
        'NaN right-hand side reports failure')

    ! 6. Status is reset by the next successful call.
    y1 = 1.0_dp
    call odeint_allroutines_checked(y1, 1, 0.0_dp, 0.5_dp, rtol, rhs_blowup, ierr)
    err = abs(y1(1) - 2.0_dp)
    call check(ierr == 0 .and. .not. odeint_has_failed() .and. err < 1.0e-8_dp, &
        'status reset after success')

    ! A non-finite initial state cannot be accepted even when the RHS is finite.
    y1 = ieee_value(1.0_dp, ieee_quiet_nan)
    call odeint_allroutines_checked(y1, 1, 0.0_dp, 0.1_dp, rtol, rhs_constant, ierr)
    call check(ierr == 1 .and. odeint_has_failed(), &
        'non-finite initial state reports failure')

    y1 = ieee_value(1.0_dp, ieee_positive_inf)
    call odeint_allroutines_checked(y1, 1, 0.0_dp, 0.1_dp, rtol, rhs_constant, ierr)
    call check(ierr /= 0, 'infinite initial state reports failure')
    y1 = 0.75_dp*huge(1.0_dp)
    call odeint_allroutines_checked(y1, 1, 0.0_dp, 1.0_dp, rtol, &
        rhs_large_constant, ierr, atol=[rtol])
    call check(ierr /= 0, 'overflowing trial state reports failure')

    ! Tolerance magnitude cannot permit a nonfinite trial endpoint.
    y1 = 0.75_dp*huge(1.0_dp)
    call odeint_allroutines_checked(y1, 1, 0.0_dp, 1.0_dp, 1.0e20_dp, &
        rhs_large_constant, ierr)
    call check(ierr /= 0, 'large relative tolerance rejects overflow')
    y1 = 0.0_dp
    call odeint_allroutines_checked(y1, 1, 0.0_dp, 1.0_dp, 1.0e20_dp, &
        rhs_nan_after_half, ierr)
    call check(ierr /= 0, 'large relative tolerance rejects NaN RHS')
    y1 = 0.0_dp
    call odeint_allroutines_checked(y1, 1, 0.0_dp, 1.0_dp, 1.0e20_dp, &
        rhs_constant, ierr)
    call check(ierr == 0 .and. abs(y1(1) - 1.0_dp) < 1.0e-14_dp, &
        'large finite relative tolerance preserves finite endpoint')
    y1 = 0.0_dp
    call odeint_allroutines(y1, 1, 0, 0.0_dp, 1.0_dp, 1.0e20_dp, &
        rhs_nan_ctx, ierr=ierr)
    call check(ierr /= 0, 'context rejects NaN RHS at large relative tolerance')

    ! Infinity is invalid even for a zero interval; finite huge remains allowed.
    y1 = 0.0_dp
    call odeint_allroutines_checked(y1, 1, 0.0_dp, 0.0_dp, rtol, &
        rhs_constant, ierr, &
        atol=[ieee_value(1.0_dp, ieee_positive_inf)])
    call check(ierr /= 0, 'zero span rejects infinite absolute tolerance')
    y1 = 0.0_dp
    call odeint_allroutines_checked(y1, 1, 0.0_dp, 1.0_dp, rtol, &
        rhs_constant, ierr, &
        atol=[ieee_value(1.0_dp, ieee_positive_inf)])
    call check(ierr /= 0, 'nonzero span rejects infinite absolute tolerance')
    y1 = 0.0_dp
    call odeint_allroutines(y1, 1, 0, 0.0_dp, 0.0_dp, rtol, &
        rhs_nan_ctx, ierr=ierr, &
        atol=[ieee_value(1.0_dp, ieee_positive_inf)])
    call check(ierr /= 0, 'context rejects infinite absolute tolerance')
    y1 = 0.0_dp
    call odeint_allroutines_checked(y1, 1, 0.0_dp, 1.0_dp, rtol, &
        rhs_constant, ierr, atol=[huge(1.0_dp)])
    call check(ierr == 0 .and. abs(y1(1) - 1.0_dp) < 1.0e-14_dp, &
        'finite huge absolute tolerance preserves analytic endpoint')

    ! 7. Context variant: same behaviour, status through ierr.
    y2 = 0.0_dp
    call odeint_allroutines(y2, 2, 0, 0.0_dp, 1.0_dp, rtol, rhs_switch_on_ctx, &
        ierr=ierr)
    call check(ierr == 1, 'context variant reports failure')
    y2 = 0.0_dp
    call odeint_allroutines(y2, 2, 0, 0.0_dp, 1.0_dp, rtol, rhs_switch_on_ctx, &
        ierr=ierr, atol=[rtol, rtol])
    err = max(abs(y2(1) - sin(1.0_dp)), abs(y2(2) - 1.0_dp))
    call check(ierr == 0 .and. err < 1.0e-8_dp, 'context variant with atol')

    if (failed) error stop 'test_odeint_abs_tol failed'
    write (*, *) 'All odeint absolute tolerance tests passed!'

contains

    subroutine check(ok, name)
        logical, intent(in) :: ok
        character(*), intent(in) :: name
        if (ok) then
            write (*, *) 'PASS: ', name
        else
            write (*, *) 'FAIL: ', name
            failed = .true.
        end if
    end subroutine check
end program test_odeint_abs_tol
