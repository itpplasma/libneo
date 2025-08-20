program test_golden_record_odeint
    !> Golden record test comparing new tableau-based implementation 
    !> against the last tagged version (v2025.07.30) with machine precision
    !>
    !> This test dynamically fetches the old implementation and compares
    !> results on multiple test problems to ensure numerical equivalence.
    
    use odeint_allroutines_sub, only: odeint_allroutines
    use odeint_golden_sub, only: odeint_allroutines_golden => odeint_allroutines
    use odeint_test_common, only: dp, linear_rhs, exponential_rhs, harmonic_rhs
    implicit none
    
    real(dp), parameter :: tol = epsilon(1.0_dp) * 200  ! 200x machine epsilon for optimized implementation comparison
    
    write(*, '(A)') '=== Golden Record Test: Tableau vs Tagged Version ==='
    write(*, '(A)') 'Comparing new implementation against v2025.07.30'
    write(*, '(A)') ''
    
    ! Run comprehensive test suite
    call test_linear_ode()
    call test_exponential_decay() 
    call test_harmonic_oscillator()
    call test_van_der_pol()
    call test_lorenz_attractor()
    call test_stiff_system()
    
    write(*, '(A)') ''
    write(*, '(A)') '✅ All golden record tests PASSED - numerical equivalence verified'
    write(*, '(A,ES10.3)') 'Maximum tolerance used: ', tol
    write(*, '(A)') ''
    
    ! Run performance benchmarks
    write(*, '(A)') '=== Performance Benchmarks ==='
    call benchmark_implementations()
    call benchmark_context_overhead()
    
contains

    subroutine test_linear_ode()
        !> Test: dy/dx = 2x, y(0) = 1, exact solution: y = x² + 1
        real(dp) :: y_new(1), y_old(1)
        real(dp), parameter :: x1 = 0.0_dp, x2 = 1.0_dp, eps = 1.0e-8_dp
        
        write(*, '(A)', advance='no') ' Testing linear ODE (dy/dx = 2x)...'
        
        y_new = [1.0_dp]
        y_old = [1.0_dp]
        
        call odeint_allroutines(y_new, 1, x1, x2, eps, linear_rhs)
        call odeint_allroutines_golden(y_old, 1, x1, x2, eps, linear_rhs)
        
        if (abs(y_new(1) - y_old(1)) > tol) then
            write(*, '(A)') ' FAILED'
            write(*, '(A,ES15.8,A,ES15.8)') '  New:', y_new(1), ' Old:', y_old(1)
            error stop 'Linear ODE golden record test failed'
        end if
        write(*, '(A)') ' PASSED'
    end subroutine test_linear_ode

    subroutine test_exponential_decay()
        !> Test: dy/dx = -2y, y(0) = 1, exact solution: y = exp(-2x)
        real(dp) :: y_new(1), y_old(1)
        real(dp), parameter :: x1 = 0.0_dp, x2 = 2.0_dp, eps = 1.0e-10_dp
        
        write(*, '(A)', advance='no') ' Testing exponential decay (dy/dx = -2y)...'
        
        y_new = [1.0_dp]
        y_old = [1.0_dp]
        
        call odeint_allroutines(y_new, 1, x1, x2, eps, exponential_rhs)
        call odeint_allroutines_golden(y_old, 1, x1, x2, eps, exponential_rhs)
        
        if (abs(y_new(1) - y_old(1)) > tol) then
            write(*, '(A)') ' FAILED'
            write(*, '(A,ES15.8,A,ES15.8)') '  New:', y_new(1), ' Old:', y_old(1)
            error stop 'Exponential decay golden record test failed'
        end if
        write(*, '(A)') ' PASSED'
    end subroutine test_exponential_decay

    subroutine test_harmonic_oscillator()
        !> Test: d²y/dt² + y = 0 → dy₁/dt = y₂, dy₂/dt = -y₁
        real(dp) :: y_new(2), y_old(2)
        real(dp), parameter :: x1 = 0.0_dp, x2 = 6.28318530717959_dp, eps = 1.0e-10_dp
        
        write(*, '(A)', advance='no') ' Testing harmonic oscillator...'
        
        y_new = [1.0_dp, 0.0_dp]  ! y(0)=1, y prime(0)=0
        y_old = [1.0_dp, 0.0_dp]
        
        call odeint_allroutines(y_new, 2, x1, x2, eps, harmonic_rhs)
        call odeint_allroutines_golden(y_old, 2, x1, x2, eps, harmonic_rhs)
        
        if (maxval(abs(y_new - y_old)) > tol) then
            write(*, '(A)') ' FAILED'
            write(*, '(A,2ES15.8)') '  New:', y_new
            write(*, '(A,2ES15.8)') '  Old:', y_old
            error stop 'Harmonic oscillator golden record test failed'
        end if
        write(*, '(A)') ' PASSED'
    end subroutine test_harmonic_oscillator

    subroutine test_van_der_pol()
        !> Test: Van der Pol oscillator with μ=1.0
        real(dp) :: y_new(2), y_old(2)
        real(dp), parameter :: x1 = 0.0_dp, x2 = 10.0_dp, eps = 1.0e-8_dp
        
        write(*, '(A)', advance='no') ' Testing Van der Pol oscillator...'
        
        y_new = [2.0_dp, 0.0_dp]
        y_old = [2.0_dp, 0.0_dp]
        
        call odeint_allroutines(y_new, 2, x1, x2, eps, van_der_pol_deriv)
        call odeint_allroutines_golden(y_old, 2, x1, x2, eps, van_der_pol_deriv)
        
        if (maxval(abs(y_new - y_old)) > tol) then
            write(*, '(A)') ' FAILED'
            write(*, '(A,2ES15.8)') '  New:', y_new
            write(*, '(A,2ES15.8)') '  Old:', y_old
            error stop 'Van der Pol golden record test failed'
        end if
        write(*, '(A)') ' PASSED'
    end subroutine test_van_der_pol

    subroutine test_lorenz_attractor()
        !> Test: Lorenz system (σ=10, ρ=28, β=8/3)
        real(dp) :: y_new(3), y_old(3)
        real(dp), parameter :: x1 = 0.0_dp, x2 = 1.0_dp, eps = 1.0e-10_dp
        
        write(*, '(A)', advance='no') ' Testing Lorenz attractor...'
        
        y_new = [1.0_dp, 1.0_dp, 1.0_dp]
        y_old = [1.0_dp, 1.0_dp, 1.0_dp]
        
        call odeint_allroutines(y_new, 3, x1, x2, eps, lorenz_deriv)
        call odeint_allroutines_golden(y_old, 3, x1, x2, eps, lorenz_deriv)
        
        if (maxval(abs(y_new - y_old)) > tol) then
            write(*, '(A)') ' FAILED'
            write(*, '(A,3ES15.8)') '  New:', y_new
            write(*, '(A,3ES15.8)') '  Old:', y_old
            error stop 'Lorenz attractor golden record test failed'
        end if
        write(*, '(A)') ' PASSED'
    end subroutine test_lorenz_attractor

    subroutine test_stiff_system()
        !> Test: Stiff system dy₁/dx = -1000y₁ + 1000y₂, dy₂/dx = y₁ - 2y₂
        real(dp) :: y_new(2), y_old(2)
        real(dp), parameter :: x1 = 0.0_dp, x2 = 0.01_dp, eps = 1.0e-6_dp
        
        write(*, '(A)', advance='no') ' Testing stiff system...'
        
        y_new = [1.0_dp, 0.0_dp]
        y_old = [1.0_dp, 0.0_dp]
        
        call odeint_allroutines(y_new, 2, x1, x2, eps, stiff_deriv)
        call odeint_allroutines_golden(y_old, 2, x1, x2, eps, stiff_deriv)
        
        if (maxval(abs(y_new - y_old)) > tol) then
            write(*, '(A)') ' FAILED'
            write(*, '(A,2ES15.8)') '  New:', y_new
            write(*, '(A,2ES15.8)') '  Old:', y_old
            error stop 'Stiff system golden record test failed'
        end if
        write(*, '(A)') ' PASSED'
    end subroutine test_stiff_system

    !> Derivative functions for test problems
    !> Unique RHS functions for golden record tests
    
    subroutine van_der_pol_deriv(x, y, dydx)
        real(dp), intent(in) :: x, y(:)
        real(dp), intent(out) :: dydx(:)
        real(dp), parameter :: mu = 1.0_dp
        
        associate(x_unused => x)
        end associate
        
        dydx(1) = y(2)
        dydx(2) = mu * (1.0_dp - y(1)**2) * y(2) - y(1)
    end subroutine van_der_pol_deriv

    subroutine lorenz_deriv(x, y, dydx)
        real(dp), intent(in) :: x, y(:)
        real(dp), intent(out) :: dydx(:)
        real(dp), parameter :: sigma = 10.0_dp, rho = 28.0_dp, beta = 8.0_dp/3.0_dp
        
        associate(x_unused => x)
        end associate
        
        dydx(1) = sigma * (y(2) - y(1))
        dydx(2) = y(1) * (rho - y(3)) - y(2)
        dydx(3) = y(1) * y(2) - beta * y(3)
    end subroutine lorenz_deriv

    subroutine stiff_deriv(x, y, dydx)
        real(dp), intent(in) :: x, y(:)
        real(dp), intent(out) :: dydx(:)
        
        associate(x_unused => x)
        end associate
        
        dydx(1) = -1000.0_dp * y(1) + 1000.0_dp * y(2)
        dydx(2) = y(1) - 2.0_dp * y(2)
    end subroutine stiff_deriv

    !> Performance benchmarking subroutines
    subroutine benchmark_implementations()
        !> Compare performance: new tableau vs old implementation
        integer, parameter :: n_trials = 1000
        real(dp) :: y_new(3), y_old(3)
        real(dp), parameter :: x1 = 0.0_dp, x2 = 10.0_dp, eps = 1.0e-8_dp
        real(dp), parameter :: max_slowdown = 1.50_dp  ! 50% threshold
        real(dp) :: start_time, end_time, time_new, time_old
        integer :: i
        
        write(*, '(A)') ' Benchmarking: New Tableau vs Old Implementation'
        write(*, '(A,I0,A)') ' Running ', n_trials, ' iterations of Lorenz attractor integration'
        
        ! Benchmark new implementation
        call cpu_time(start_time)
        do i = 1, n_trials
            y_new = [1.0_dp, 1.0_dp, 1.0_dp]
            call odeint_allroutines(y_new, 3, x1, x2, eps, lorenz_deriv)
        end do
        call cpu_time(end_time)
        time_new = end_time - start_time
        
        ! Benchmark old implementation
        call cpu_time(start_time)
        do i = 1, n_trials
            y_old = [1.0_dp, 1.0_dp, 1.0_dp]
            call odeint_allroutines_golden(y_old, 3, x1, x2, eps, lorenz_deriv)
        end do
        call cpu_time(end_time)
        time_old = end_time - start_time
        
        write(*, '(A,F8.4,A)') ' New tableau implementation: ', time_new, ' seconds'
        write(*, '(A,F8.4,A)') ' Old manual implementation: ', time_old, ' seconds'
        write(*, '(A,F6.2,A)') ' Performance ratio (old/new): ', time_old/time_new, 'x'
        if (time_new < time_old) then
            write(*, '(A,F5.1,A)') ' 🚀 New implementation is ', (time_old/time_new - 1.0_dp)*100, '% faster'
        else
            write(*, '(A,F5.1,A)') ' 📉 New implementation is ', (time_new/time_old - 1.0_dp)*100, '% slower'
        end if
        
        ! Performance threshold check: fail if >50% slower
        if (time_new > time_old * max_slowdown) then
            write(*, '(A)') ''
            write(*, '(A)') '❌ PERFORMANCE REGRESSION: New implementation exceeds 50% slowdown threshold'
            write(*, '(A,F5.1,A)') ' Maximum allowed slowdown: 50.0%'
            write(*, '(A,F5.1,A)') ' Actual slowdown: ', (time_new/time_old - 1.0_dp)*100, '%'
            error stop 'Performance benchmark failed - exceeds 50% slowdown threshold'
        else
            write(*, '(A)') ' ✅ Performance within acceptable threshold (<50% slowdown)'
        end if
        write(*, '(A)') ''
    end subroutine benchmark_implementations

    subroutine benchmark_context_overhead()
        !> Compare performance: context vs no-context versions
        integer, parameter :: n_trials = 1000
        real(dp) :: y_no_ctx(3), y_ctx(3)
        real(dp), parameter :: x1 = 0.0_dp, x2 = 5.0_dp, eps = 1.0e-8_dp
        real(dp) :: start_time, end_time, time_no_ctx, time_ctx
        integer :: i
        type :: dummy_context_t
            real(dp) :: dummy_value = 1.0_dp
        end type dummy_context_t
        type(dummy_context_t) :: ctx
        
        write(*, '(A)') ' Benchmarking: Context vs No-Context Overhead'
        write(*, '(A,I0,A)') ' Running ', n_trials, ' iterations of Van der Pol integration'
        
        ! Benchmark no-context version
        call cpu_time(start_time)
        do i = 1, n_trials
            y_no_ctx = [2.0_dp, 0.0_dp, 0.0_dp] ! 3rd element unused
            call odeint_allroutines(y_no_ctx(1:2), 2, x1, x2, eps, van_der_pol_deriv)
        end do
        call cpu_time(end_time)
        time_no_ctx = end_time - start_time
        
        ! Benchmark context version  
        call cpu_time(start_time)
        do i = 1, n_trials
            y_ctx = [2.0_dp, 0.0_dp, 0.0_dp] ! 3rd element unused
            call odeint_allroutines(y_ctx(1:2), 2, ctx, x1, x2, eps, van_der_pol_context_deriv)
        end do
        call cpu_time(end_time)
        time_ctx = end_time - start_time
        
        write(*, '(A,F8.4,A)') ' No-context version: ', time_no_ctx, ' seconds'
        write(*, '(A,F8.4,A)') ' Context version: ', time_ctx, ' seconds'
        write(*, '(A,F6.2,A)') ' Overhead ratio (ctx/no-ctx): ', time_ctx/time_no_ctx, 'x'
        write(*, '(A,F5.1,A)') ' Context overhead: ', (time_ctx/time_no_ctx - 1.0_dp)*100, '%'
        if (abs(time_ctx/time_no_ctx - 1.0_dp) < 0.05_dp) then
            write(*, '(A)') ' ✅ Minimal overhead - context parameter is efficiently handled'
        else if (time_ctx/time_no_ctx > 1.1_dp) then
            write(*, '(A)') ' ⚠️  Significant overhead detected'
        end if
        write(*, '(A)') ''
    end subroutine benchmark_context_overhead
    
    !> Context-enabled derivative for overhead testing
    subroutine van_der_pol_context_deriv(x, y, dydx, context)
        real(dp), intent(in) :: x, y(:)
        real(dp), intent(out) :: dydx(:)
        class(*), intent(in) :: context
        real(dp), parameter :: mu = 1.0_dp
        
        associate(x_unused => x, context_unused => context)
        end associate
        
        dydx(1) = y(2)
        dydx(2) = mu * (1.0_dp - y(1)**2) * y(2) - y(1)
    end subroutine van_der_pol_context_deriv

end program test_golden_record_odeint