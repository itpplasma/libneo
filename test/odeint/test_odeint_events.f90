program test_odeint_events
   use libneo_kinds, only: dp
   use odeint_allroutines_sub, only: odeint_allroutines, ode_event_t
   implicit none

   type :: threshold_t
      real(dp) :: value
   end type threshold_t

   logical :: test_failed

   test_failed = .false.

   call test_terminal_event(test_failed)
   call test_direction_filter(test_failed)
   call test_non_terminal_event(test_failed)
   call test_context_event(test_failed)

   if (test_failed) then
      write (*, *) 'Some odeint event tests failed!'
      stop 1
   end if
   write (*, *) 'All odeint event tests passed!'

contains

   subroutine test_terminal_event(test_failed)
      logical, intent(inout) :: test_failed

      real(dp), parameter :: x1 = 0.0_dp, x2 = 1.0_dp, eps = 1.0e-8_dp
      real(dp) :: y(1)
      type(ode_event_t) :: ev(1)

      ev(1)%f => event_y_minus_half
      ev(1)%direction = 0
      ev(1)%terminal = .true.

      y(1) = 0.0_dp
      call odeint_allroutines(y, 1, x1, x2, eps, rhs_unit, events=ev)

      if (.not. ev(1)%found) then
         write (*, *) 'ERROR: terminal event not detected'
         test_failed = .true.
      end if
      if (abs(ev(1)%x - 0.5_dp) > 1.0e-6_dp) then
         write (*, *) 'ERROR: terminal event time incorrect', ev(1)%x
         test_failed = .true.
      end if
      if (abs(y(1) - 0.5_dp) > 1.0e-6_dp) then
         write (*, *) 'ERROR: terminal event state incorrect', y(1)
         test_failed = .true.
      end if
   end subroutine test_terminal_event

   subroutine test_direction_filter(test_failed)
      logical, intent(inout) :: test_failed

      real(dp), parameter :: x1 = 0.0_dp, x2 = 1.0_dp, eps = 1.0e-8_dp
      real(dp) :: y(1)
      type(ode_event_t) :: ev(1)

      ev(1)%f => event_y_minus_half
      ev(1)%direction = -1
      ev(1)%terminal = .true.

      y(1) = 0.0_dp
      call odeint_allroutines(y, 1, x1, x2, eps, rhs_unit, events=ev)

      if (ev(1)%found) then
         write (*, *) 'ERROR: direction filter should block event'
         test_failed = .true.
      end if
      if (abs(y(1) - 1.0_dp) > 1.0e-6_dp) then
         write (*, *) 'ERROR: direction filter end state incorrect', y(1)
         test_failed = .true.
      end if
   end subroutine test_direction_filter

   subroutine test_non_terminal_event(test_failed)
      logical, intent(inout) :: test_failed

      real(dp), parameter :: x1 = 0.0_dp, x2 = 1.0_dp, eps = 1.0e-8_dp
      real(dp) :: y(1)
      type(ode_event_t) :: ev(1)

      ev(1)%f => event_y_minus_half
      ev(1)%direction = 0
      ev(1)%terminal = .false.

      y(1) = 0.0_dp
      call odeint_allroutines(y, 1, x1, x2, eps, rhs_unit, events=ev)

      if (.not. ev(1)%found) then
         write (*, *) 'ERROR: non-terminal event not detected'
         test_failed = .true.
      end if
      if (abs(y(1) - 1.0_dp) > 1.0e-6_dp) then
         write (*, *) 'ERROR: non-terminal event should continue', y(1)
         test_failed = .true.
      end if
   end subroutine test_non_terminal_event

   subroutine test_context_event(test_failed)
      logical, intent(inout) :: test_failed

      real(dp), parameter :: x1 = 0.0_dp, x2 = 1.0_dp, eps = 1.0e-8_dp
      real(dp) :: y(1)
      type(ode_event_t) :: ev(1)
      type(threshold_t) :: ctx

      ctx%value = 0.25_dp
      ev(1)%f => event_with_context
      ev(1)%direction = 0
      ev(1)%terminal = .true.

      y(1) = 0.0_dp
      call odeint_allroutines(y, 1, ctx, x1, x2, eps, rhs_unit_ctx, events=ev)

      if (.not. ev(1)%found) then
         write (*, *) 'ERROR: context event not detected'
         test_failed = .true.
      end if
      if (abs(ev(1)%x - 0.25_dp) > 1.0e-6_dp) then
         write (*, *) 'ERROR: context event time incorrect', ev(1)%x
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

   subroutine event_y_minus_half(x, y, g, context)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: g
      class(*), intent(in), optional :: context

      associate (x_unused => x)
      end associate
      associate (context_unused => context)
      end associate

      g = y(1) - 0.5_dp
   end subroutine event_y_minus_half

   subroutine event_with_context(x, y, g, context)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: g
      class(*), intent(in), optional :: context

      associate (x_unused => x)
      end associate

      select type (context)
      type is (threshold_t)
         g = y(1) - context%value
      class default
         g = y(1)
      end select
   end subroutine event_with_context

end program test_odeint_events
