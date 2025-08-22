program test_odeint_allroutines_context
    use odeint_test_common
    implicit none

    logical :: test_failed = .false.

    call run_all_serial_tests(test_failed)

    if (test_failed) then
        write(*,*) 'Some odeint_allroutines tests failed!'
        error stop 'Serial test failed'
    else
        write(*,*) 'All odeint_allroutines tests passed!'
    end if

contains

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