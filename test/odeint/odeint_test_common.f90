module odeint_test_common
    use libneo_kinds, only: dp
    use odeint_allroutines_sub, only: odeint_allroutines
    implicit none

    type :: ode_params_t
        real(dp) :: omega = 1.0_dp
        real(dp) :: damping = 0.1_dp
        real(dp) :: amplitude = 2.0_dp
    end type ode_params_t

    real(dp), parameter :: TOL_LOOSE = 1e-4_dp
    real(dp), parameter :: TOL_TIGHT = 1e-5_dp
    real(dp), parameter :: TOL_VERY_TIGHT = 1e-6_dp

contains

    ! Test runner that can be used by both serial and parallel tests
    subroutine run_odeint_test(test_name, test_type, thread_id, iteration, test_failed)
        character(len=*), intent(in) :: test_name
        integer, intent(in) :: test_type, thread_id, iteration
        logical, intent(inout) :: test_failed
        
        select case (test_type)
        case (1)
            call test_linear_ode_impl(thread_id, iteration, test_failed)
        case (2)
            call test_exponential_decay_impl(thread_id, iteration, test_failed)
        case (3)
            call test_harmonic_oscillator_impl(thread_id, iteration, test_failed)
        case (4)
            call test_multi_variable_system_impl(thread_id, iteration, test_failed)
        case (5)
            call test_context_parameters_impl(thread_id, iteration, test_failed)
        case (6)
            call test_initial_stepsize_impl(thread_id, iteration, test_failed)
        case (7)
            call test_backward_integration_impl(thread_id, iteration, test_failed)
        case default
            !$omp critical
            write(*,*) 'ERROR: Unknown test type', test_type
            test_failed = .true.
            !$omp end critical
        end select
    end subroutine run_odeint_test

    subroutine test_linear_ode_impl(thread_id, iteration, test_failed)
        integer, intent(in) :: thread_id, iteration
        logical, intent(inout) :: test_failed
        
        real(dp), parameter :: x1 = 0.0_dp, x2 = 1.0_dp, eps = 1.0e-8_dp
        integer, parameter :: nvar = 1
        real(dp) :: y(nvar)
        
        y(1) = 0.0_dp
        call odeint_allroutines(y, nvar, x1, x2, eps, linear_rhs)
        
        if (abs(y(1) - 1.0_dp) > TOL_TIGHT) then
            !$omp critical
            write(*,*) 'ERROR: Linear ODE failed on thread', thread_id, 'iteration', iteration
            write(*,*) 'Expected 1.0, got', y(1)
            test_failed = .true.
            !$omp end critical
        end if
    end subroutine test_linear_ode_impl

    subroutine test_exponential_decay_impl(thread_id, iteration, test_failed)
        integer, intent(in) :: thread_id, iteration
        logical, intent(inout) :: test_failed
        
        real(dp), parameter :: x1 = 0.0_dp, x2 = 2.0_dp, eps = 1.0e-8_dp
        integer, parameter :: nvar = 1
        real(dp) :: y(nvar)
        real(dp) :: exact_solution
        
        y(1) = 1.0_dp
        call odeint_allroutines(y, nvar, x1, x2, eps, exponential_rhs)
        
        exact_solution = exp(-2.0_dp * x2)
        if (abs(y(1) - exact_solution) > TOL_TIGHT) then
            !$omp critical
            write(*,*) 'ERROR: Exponential decay failed on thread', thread_id, 'iteration', iteration
            write(*,*) 'Expected', exact_solution, 'got', y(1)
            test_failed = .true.
            !$omp end critical
        end if
    end subroutine test_exponential_decay_impl

    subroutine test_harmonic_oscillator_impl(thread_id, iteration, test_failed)
        integer, intent(in) :: thread_id, iteration
        logical, intent(inout) :: test_failed
        
        real(dp), parameter :: x1 = 0.0_dp, x2 = 1.0_dp, eps = 1.0e-8_dp
        integer, parameter :: nvar = 2
        real(dp) :: y(nvar)
        real(dp) :: expected_x, expected_v
        
        y(1) = 1.0_dp  ! position
        y(2) = 0.0_dp  ! velocity
        call odeint_allroutines(y, nvar, x1, x2, eps, harmonic_rhs)
        
        expected_x = cos(x2)
        expected_v = -sin(x2)
        
        if (abs(y(1) - expected_x) > TOL_LOOSE .or. abs(y(2) - expected_v) > TOL_LOOSE) then
            !$omp critical
            write(*,*) 'ERROR: Harmonic oscillator failed on thread', thread_id, 'iteration', iteration
            write(*,*) 'Position: expected', expected_x, 'got', y(1)
            write(*,*) 'Velocity: expected', expected_v, 'got', y(2)
            test_failed = .true.
            !$omp end critical
        end if
    end subroutine test_harmonic_oscillator_impl

    subroutine test_multi_variable_system_impl(thread_id, iteration, test_failed)
        integer, intent(in) :: thread_id, iteration
        logical, intent(inout) :: test_failed
        
        real(dp), parameter :: x1 = 0.0_dp, x2 = 0.5_dp, eps = 1.0e-8_dp
        integer, parameter :: nvar = 3
        real(dp) :: y(nvar)
        
        y(1) = 1.0_dp
        y(2) = 0.0_dp
        y(3) = -1.0_dp
        call odeint_allroutines(y, nvar, x1, x2, eps, coupled_system_rhs)
        
        if (any(abs(y) > 10.0_dp)) then
            !$omp critical
            write(*,*) 'ERROR: Coupled system solution exploded on thread', thread_id, 'iteration', iteration
            write(*,*) 'Solution:', y
            test_failed = .true.
            !$omp end critical
        end if
    end subroutine test_multi_variable_system_impl

    subroutine test_context_parameters_impl(thread_id, iteration, test_failed)
        integer, intent(in) :: thread_id, iteration
        logical, intent(inout) :: test_failed
        
        real(dp), parameter :: x1 = 0.0_dp, x2 = 1.0_dp, eps = 1.0e-8_dp
        integer, parameter :: nvar = 2
        real(dp) :: y(nvar)
        type(ode_params_t) :: params
        
        ! Vary parameters slightly for each iteration to test different scenarios
        params%omega = 2.0_dp + 0.1_dp * mod(iteration, 10)
        params%damping = 0.2_dp + 0.05_dp * mod(iteration, 5)
        
        y(1) = 1.0_dp  ! position
        y(2) = 0.0_dp  ! velocity
        call odeint_allroutines(y, nvar, params, x1, x2, eps, damped_oscillator_rhs)
        
        if (abs(y(1)) > 1.0_dp) then
            !$omp critical
            write(*,*) 'ERROR: Damped oscillator amplitude too large on thread', thread_id, 'iteration', iteration
            write(*,*) 'Amplitude:', y(1)
            test_failed = .true.
            !$omp end critical
        end if
    end subroutine test_context_parameters_impl

    subroutine test_initial_stepsize_impl(thread_id, iteration, test_failed)
        integer, intent(in) :: thread_id, iteration
        logical, intent(inout) :: test_failed
        
        real(dp), parameter :: x1 = 0.0_dp, x2 = 0.1_dp, eps = 1.0e-8_dp
        real(dp), parameter :: h_init = 0.001_dp
        integer, parameter :: nvar = 1
        real(dp) :: y(nvar)
        
        y(1) = 1.0_dp
        call odeint_allroutines(y, nvar, x1, x2, eps, exponential_rhs, h_init)
        
        if (abs(y(1) - exp(-2.0_dp * x2)) > TOL_VERY_TIGHT) then
            !$omp critical
            write(*,*) 'ERROR: Initial stepsize test failed on thread', thread_id, 'iteration', iteration
            test_failed = .true.
            !$omp end critical
        end if
    end subroutine test_initial_stepsize_impl

    subroutine test_backward_integration_impl(thread_id, iteration, test_failed)
        integer, intent(in) :: thread_id, iteration
        logical, intent(inout) :: test_failed
        
        real(dp), parameter :: x1 = 1.0_dp, x2 = 0.0_dp, eps = 1.0e-8_dp
        integer, parameter :: nvar = 1
        real(dp) :: y(nvar)
        
        y(1) = exp(-2.0_dp)
        call odeint_allroutines(y, nvar, x1, x2, eps, exponential_rhs)
        
        if (abs(y(1) - 1.0_dp) > TOL_TIGHT) then
            !$omp critical
            write(*,*) 'ERROR: Backward integration failed on thread', thread_id, 'iteration', iteration
            write(*,*) 'Expected 1.0, got', y(1)
            test_failed = .true.
            !$omp end critical
        end if
    end subroutine test_backward_integration_impl

    ! Common RHS functions
    subroutine linear_rhs(x, y, dydx)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydx(:)
        
        associate(y_unused => y)
        end associate
        
        dydx(1) = 2.0_dp * x
    end subroutine linear_rhs

    subroutine exponential_rhs(x, y, dydx)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydx(:)
        
        associate(x_unused => x)
        end associate
        
        dydx(1) = -2.0_dp * y(1)
    end subroutine exponential_rhs

    subroutine harmonic_rhs(x, y, dydx)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydx(:)
        
        associate(x_unused => x)
        end associate
        
        dydx(1) = y(2)      ! dx/dt = v
        dydx(2) = -y(1)     ! dv/dt = -x
    end subroutine harmonic_rhs

    subroutine coupled_system_rhs(x, y, dydx)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydx(:)
        
        associate(x_unused => x)
        end associate
        
        dydx(1) = -0.1_dp * y(1) + 0.05_dp * y(2)
        dydx(2) = 0.05_dp * y(1) - 0.1_dp * y(2) + 0.02_dp * y(3)
        dydx(3) = 0.02_dp * y(2) - 0.05_dp * y(3)
    end subroutine coupled_system_rhs

    subroutine damped_oscillator_rhs(x, y, dydx, context)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydx(:)
        class(*), intent(in) :: context
        
        associate(x_unused => x)
        end associate
        
        select type (context)
        type is (ode_params_t)
            dydx(1) = y(2)
            dydx(2) = -context%omega**2 * y(1) - 2.0_dp * context%damping * context%omega * y(2)
        end select
    end subroutine damped_oscillator_rhs

end module odeint_test_common