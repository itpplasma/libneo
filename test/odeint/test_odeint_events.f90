program test_odeint_events
   use libneo_kinds, only: dp
   use odeint_allroutines_sub, only: odeint_allroutines, ode_event_t
   implicit none

   type :: threshold_t
      real(dp) :: value
   end type threshold_t

   logical :: test_failed

   test_failed = .false.

   call test_event_trigger(test_failed)
   call test_event_approx_start(test_failed)
   call test_direction_filter(test_failed)
   call test_context_event(test_failed)

   if (test_failed) then
      write (*, *) 'Some odeint event tests failed!'
      stop 1
   end if
   write (*, *) 'All odeint event tests passed!'

contains

   subroutine test_event_trigger(test_failed)
      logical, intent(inout) :: test_failed

      real(dp), parameter :: x1 = 0.0_dp, x2 = 1.0_dp, eps = 1.0e-8_dp
      real(dp) :: y(1)
      type(ode_event_t) :: ev(1)

      ev(1)%condition => event_y_minus_half
      ev(1)%direction = 0
      y(1) = 0.0_dp
      call odeint_allroutines(y, 1, x1, x2, eps, rhs_unit, events=ev)

      if (.not. ev(1)%triggered) then
         write (*, *) 'ERROR: event not detected'
         test_failed = .true.
      end if
      if (abs(ev(1)%x_event - 0.5_dp) > 1.0e-6_dp) then
         write (*, *) 'ERROR: event time incorrect', ev(1)%x_event
         test_failed = .true.
      end if
      if (abs(y(1) - 0.5_dp) > 1.0e-6_dp) then
         write (*, *) 'ERROR: event state incorrect', y(1)
         test_failed = .true.
      end if
   end subroutine test_event_trigger

   subroutine test_event_approx_start(test_failed)
      logical, intent(inout) :: test_failed

      real(dp), parameter :: x1 = 0.0_dp, x2 = 1.0_dp, eps = 1.0e-8_dp
      real(dp) :: y(1)
      type(ode_event_t) :: ev(1)

      ev(1)%condition => event_y
      ev(1)%direction = 0
      y(1) = 1.0e-10_dp
      call odeint_allroutines(y, 1, x1, x2, eps, rhs_unit, events=ev)

      if (.not. ev(1)%triggered) then
         write (*, *) 'ERROR: approx-start event not detected'
         test_failed = .true.
      end if
      if (abs(ev(1)%x_event - x1) > 1.0e-12_dp) then
         write (*, *) 'ERROR: approx-start event time incorrect', ev(1)%x_event
         test_failed = .true.
      end if
      if (abs(y(1) - 1.0e-10_dp) > 1.0e-8_dp) then
         write (*, *) 'ERROR: approx-start event state incorrect', y(1)
         test_failed = .true.
      end if
   end subroutine test_event_approx_start

   subroutine test_direction_filter(test_failed)
      logical, intent(inout) :: test_failed

      real(dp), parameter :: x1 = 0.0_dp, x2 = 1.0_dp, eps = 1.0e-8_dp
      real(dp) :: y(1)
      type(ode_event_t) :: ev(1)

      ev(1)%condition => event_y_minus_half
      ev(1)%direction = -1
      y(1) = 0.0_dp
      call odeint_allroutines(y, 1, x1, x2, eps, rhs_unit, events=ev)

      if (ev(1)%triggered) then
         write (*, *) 'ERROR: direction filter should block event'
         test_failed = .true.
      end if
      if (abs(y(1) - 1.0_dp) > 1.0e-6_dp) then
         write (*, *) 'ERROR: direction filter end state incorrect', y(1)
         test_failed = .true.
      end if
   end subroutine test_direction_filter

   subroutine test_context_event(test_failed)
      logical, intent(inout) :: test_failed

      real(dp), parameter :: x1 = 0.0_dp, x2 = 1.0_dp, eps = 1.0e-8_dp
      real(dp) :: y(1)
      type(ode_event_t) :: ev(1)
      type(threshold_t) :: ctx

      ctx%value = 0.25_dp
      ev(1)%condition => event_with_context
      ev(1)%direction = 0
      y(1) = 0.0_dp
      call odeint_allroutines(y, 1, ctx, x1, x2, eps, rhs_unit_ctx, events=ev)

      if (.not. ev(1)%triggered) then
         write (*, *) 'ERROR: context event not detected'
         test_failed = .true.
      end if
      if (abs(ev(1)%x_event - 0.25_dp) > 1.0e-6_dp) then
         write (*, *) 'ERROR: context event time incorrect', ev(1)%x_event
         test_failed = .true.
      end if
   end subroutine test_context_event

   subroutine rhs_unit(x, y, dydx)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: dydx(:)

      associate (x_unused => x)
      end associate
      associate (y_unused => y)
      end associate

      dydx(1) = 1.0_dp
   end subroutine rhs_unit

   subroutine rhs_unit_ctx(x, y, dydx, context)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: dydx(:)
      class(*), intent(in) :: context

      associate (x_unused => x)
      end associate
      associate (y_unused => y)
      end associate

      select type (context)
      type is (threshold_t)
         dydx(1) = 1.0_dp
      end select
   end subroutine rhs_unit_ctx

   real(dp) function event_y_minus_half(x, y, context)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: y(:)
      class(*), intent(in), optional :: context

      associate (x_unused => x)
      end associate
      associate (context_unused => context)
      end associate

      event_y_minus_half = y(1) - 0.5_dp
   end function event_y_minus_half

   real(dp) function event_y(x, y, context)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: y(:)
      class(*), intent(in), optional :: context

      associate (x_unused => x)
      end associate
      associate (context_unused => context)
      end associate

      event_y = y(1)
   end function event_y

   real(dp) function event_with_context(x, y, context)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: y(:)
      class(*), intent(in), optional :: context

      associate (x_unused => x)
      end associate

      select type (context)
      type is (threshold_t)
         event_with_context = y(1) - context%value
      class default
         event_with_context = y(1)
      end select
   end function event_with_context

end program test_odeint_events
