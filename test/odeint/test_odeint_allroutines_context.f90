program test_odeint_allroutines_context
    use odeint_test_common
    use odeint_allroutines_sub, only: odeint_allroutines, odeint_has_failed
    implicit none

    logical :: test_failed = .false.

    call run_all_serial_tests(test_failed)
    call test_step_limit_status(test_failed)

    if (test_failed) then
        write(*,*) 'Some odeint_allroutines tests failed!'
        error stop 'Serial test failed'
    else
        write(*,*) 'All odeint_allroutines tests passed!'
    end if

contains

    subroutine test_step_limit_status(test_failed)
        logical, intent(inout) :: test_failed
        real(dp) :: y(2)
        type(ode_params_t) :: params
        integer :: ierr

        write(*,*) 'Testing catchable step-limit failure'
        y = 0.0_dp
        call odeint_allroutines(y, 2, params, 0.0_dp, 1.0_dp, 1.0e-12_dp, &
                                damped_oscillator_rhs, initial_stepsize=1.0e-3_dp, &
                                step_limit=1, ierr=ierr)
        if (ierr /= 1 .or. .not. odeint_has_failed()) then
            write(*,*) 'ERROR: exhausted step limit did not set failure status'
            test_failed = .true.
        else
            write(*,*) 'Step-limit status test passed'
        end if
        y = 0.0_dp
        call odeint_allroutines(y, 2, params, 0.0_dp, 1.0_dp, 1.0e-12_dp, &
                                damped_oscillator_rhs, ierr=ierr)
        if (ierr /= 0 .or. odeint_has_failed()) then
            write(*,*) 'ERROR: failure status leaked into the next solve'
            test_failed = .true.
        end if
    end subroutine test_step_limit_status

    subroutine run_all_serial_tests(test_failed)
        logical, intent(inout) :: test_failed
        
        write(*,*) 'Testing linear ODE: dy/dx = 2x'
        call run_odeint_test('Linear ODE', 1, 0, 0, test_failed)
        if (.not. test_failed) write(*,*) 'Linear ODE test passed'
        
        write(*,*) 'Testing exponential decay: dy/dx = -2y'
        call run_odeint_test('Exponential decay', 2, 0, 0, test_failed)
        if (.not. test_failed) write(*,*) 'Exponential decay test passed'
        
        write(*,*) 'Testing harmonic oscillator: d²x/dt² = -x'
        call run_odeint_test('Harmonic oscillator', 3, 0, 0, test_failed)
        if (.not. test_failed) write(*,*) 'Harmonic oscillator test passed'
        
        write(*,*) 'Testing 3-variable coupled system'
        call run_odeint_test('Multi-variable system', 4, 0, 0, test_failed)
        if (.not. test_failed) write(*,*) 'Multi-variable system test passed'
        
        write(*,*) 'Testing parameterized damped oscillator with context'
        call run_odeint_test('Context parameters', 5, 0, 0, test_failed)
        if (.not. test_failed) write(*,*) 'Context parameter test passed'
        
        write(*,*) 'Testing initial stepsize parameter'
        call run_odeint_test('Initial stepsize', 6, 0, 0, test_failed)
        if (.not. test_failed) write(*,*) 'Initial stepsize test passed'
        
        write(*,*) 'Testing backward integration'
        call run_odeint_test('Backward integration', 7, 0, 0, test_failed)
        if (.not. test_failed) write(*,*) 'Backward integration test passed'
    end subroutine run_all_serial_tests

end program test_odeint_allroutines_context
