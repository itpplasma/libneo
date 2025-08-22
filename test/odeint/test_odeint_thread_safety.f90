program test_odeint_thread_safety
    use omp_lib
    use odeint_test_common
    implicit none

    integer, parameter :: NUM_THREADS = 8
    integer, parameter :: NUM_ITERATIONS = 16
    logical :: all_tests_passed = .true.

    write(*,*) 'Testing thread safety of odeint_allroutines with', NUM_THREADS, 'threads'
    write(*,*) 'Running', NUM_ITERATIONS, 'iterations per thread'

    call omp_set_num_threads(NUM_THREADS)

    ! Test all scenarios in parallel
    call test_scenario_parallel('Linear ODE', 1)
    call test_scenario_parallel('Exponential decay', 2)
    call test_scenario_parallel('Harmonic oscillator', 3)
    call test_scenario_parallel('Context parameters', 4)
    call test_scenario_parallel('Mixed scenarios', 0)

    if (all_tests_passed) then
        write(*,*) 'All thread-safety tests passed!'
    else
        write(*,*) 'Some thread-safety tests failed!'
        error stop 'Thread-safety test failed'
    end if

contains

    subroutine test_scenario_parallel(scenario_name, test_type)
        character(len=*), intent(in) :: scenario_name
        integer, intent(in) :: test_type  ! 0=mixed, 1=linear, 2=exp, 3=harmonic, 4=context
        
        integer :: i, thread_id
        logical :: test_failed = .false.
        
        write(*,*) 'Testing', trim(scenario_name), 'in parallel...'
        
        !$omp parallel private(i, thread_id) shared(test_failed)
        thread_id = omp_get_thread_num()
        
        !$omp do
        do i = 1, NUM_ITERATIONS
            call run_test_by_type(test_type, thread_id, i, test_failed)
        end do
        !$omp end do
        !$omp end parallel
        
        if (test_failed) then
            all_tests_passed = .false.
            write(*,*) trim(scenario_name), 'parallel test failed!'
        else
            write(*,*) trim(scenario_name), 'parallel test passed'
        end if
    end subroutine test_scenario_parallel

    subroutine run_test_by_type(test_type, thread_id, iteration, test_failed)
        integer, intent(in) :: test_type, thread_id, iteration
        logical, intent(inout) :: test_failed
        
        if (test_type == 0) then
            ! Mixed scenarios - choose test based on iteration
            select case (mod(iteration, 4))
            case (0)
                call run_odeint_test('Linear ODE', 1, thread_id, iteration, test_failed)
            case (1)
                call run_odeint_test('Exponential decay', 2, thread_id, iteration, test_failed)
            case (2)
                call run_odeint_test('Harmonic oscillator', 3, thread_id, iteration, test_failed)
            case (3)
                call run_odeint_test('Context parameters', 5, thread_id, iteration, test_failed)
            end select
        else
            call run_odeint_test('Test', test_type, thread_id, iteration, test_failed)
        end if
    end subroutine run_test_by_type

end program test_odeint_thread_safety