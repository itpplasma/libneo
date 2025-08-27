program test_batch_interpolate
    use libneo_kinds, only : dp
    use math_constants
    use libneo_util
    
    implicit none
    
    real(dp), parameter :: TOL = 1.0d-12  ! Must match individual splines exactly
    real(dp), parameter :: X_MIN = 1.23d0, X_MAX = TWOPI + 1.23d0
    
    ! Test 1D batch splines
    call test_batch_spline_1d_construction()
    call test_batch_spline_1d_evaluation()
    call test_batch_spline_1d_derivatives()
    call test_batch_spline_1d_single_extraction()
    call test_batch_spline_1d_periodic()
    call test_batch_spline_1d_memory_layout()
    call bench_batch_vs_individual_1d()
    
    ! Test 2D batch splines  
    call test_batch_spline_2d_construction()
    call test_batch_spline_2d_evaluation()
    call test_batch_spline_2d_derivatives()
    call test_batch_spline_2d_mixed_periodic()
    
    ! Test 3D batch splines
    call test_batch_spline_3d_construction()
    call test_batch_spline_3d_evaluation()
    call test_batch_spline_3d_derivatives()
    call test_batch_spline_3d_field_components()
    
    print *, "All batch spline tests passed!"
    
contains

    subroutine test_batch_spline_1d_construction()
        use interpolate
        
        integer, parameter :: N_POINTS = 50
        integer, parameter :: N_QUANTITIES = 4
        integer, parameter :: ORDER = 5
        logical, parameter :: PERIODIC = .false.
        
        type(BatchSplineData1D) :: batch_spl
        type(SplineData1D) :: single_spls(N_QUANTITIES)
        
        real(dp) :: x(N_POINTS)
        real(dp) :: y_batch(N_POINTS, N_QUANTITIES)
        real(dp) :: y_single(N_POINTS)
        
        integer :: iq, i, k
        
        print *, "Testing 1D batch spline construction..."
        
        ! Generate test data
        call linspace(X_MIN, X_MAX, N_POINTS, x)
        
        do iq = 1, N_QUANTITIES
            do i = 1, N_POINTS
                y_batch(i, iq) = cos(x(i) + real(iq-1, dp)*0.5d0)
            end do
        end do
        
        ! Construct batch spline
        call construct_batch_splines_1d(X_MIN, X_MAX, y_batch, ORDER, PERIODIC, batch_spl)
        
        ! Construct individual splines for comparison
        do iq = 1, N_QUANTITIES
            y_single(:) = y_batch(:, iq)
            call construct_splines_1d(X_MIN, X_MAX, y_single, ORDER, PERIODIC, single_spls(iq))
        end do
        
        ! Verify metadata
        if (batch_spl%num_quantities /= N_QUANTITIES) error stop "Wrong num_quantities"
        if (batch_spl%num_points /= N_POINTS) error stop "Wrong num_points"
        if (batch_spl%order /= ORDER) error stop "Wrong order"
        if (batch_spl%periodic .neqv. PERIODIC) error stop "Wrong periodic flag"
        if (abs(batch_spl%x_min - X_MIN) > TOL) error stop "Wrong x_min"
        if (abs(batch_spl%h_step - single_spls(1)%h_step) > TOL) error stop "Wrong h_step"
        
        ! Verify coefficients match individual splines
        do iq = 1, N_QUANTITIES
            do i = 1, N_POINTS
                do k = 0, ORDER
                    if (abs(batch_spl%coeff(iq, k, i) - single_spls(iq)%coeff(k, i)) > TOL) then
                        print *, "Coefficient mismatch at k=", k, " i=", i, " iq=", iq
                        print *, "Batch: ", batch_spl%coeff(iq, k, i)
                        print *, "Single: ", single_spls(iq)%coeff(k, i)
                        error stop "Coefficient mismatch"
                    end if
                end do
            end do
        end do
        
        ! Clean up
        call destroy_batch_splines_1d(batch_spl)
        do iq = 1, N_QUANTITIES
            call destroy_splines_1d(single_spls(iq))
        end do
        
        print *, "  PASSED: 1D batch construction matches individual splines"
        
    end subroutine test_batch_spline_1d_construction
    
    
    subroutine test_batch_spline_1d_evaluation()
        use interpolate
        
        integer, parameter :: N_POINTS = 50
        integer, parameter :: N_QUANTITIES = 3
        integer, parameter :: ORDER = 5
        logical, parameter :: PERIODIC = .false.
        
        type(BatchSplineData1D) :: batch_spl
        type(SplineData1D) :: single_spls(N_QUANTITIES)
        
        real(dp) :: x(N_POINTS)
        real(dp) :: y_batch(N_POINTS, N_QUANTITIES)
        real(dp) :: y_single(N_POINTS)
        
        real(dp) :: x_eval
        real(dp) :: y_batch_result(N_QUANTITIES)
        real(dp) :: y_single_result
        
        integer :: iq, i, n_test
        
        print *, "Testing 1D batch spline evaluation..."
        
        ! Generate test data
        call linspace(X_MIN, X_MAX, N_POINTS, x)
        
        do iq = 1, N_QUANTITIES
            do i = 1, N_POINTS
                y_batch(i, iq) = sin(x(i)*real(iq, dp))
            end do
        end do
        
        ! Construct splines
        call construct_batch_splines_1d(X_MIN, X_MAX, y_batch, ORDER, PERIODIC, batch_spl)
        
        do iq = 1, N_QUANTITIES
            y_single(:) = y_batch(:, iq)
            call construct_splines_1d(X_MIN, X_MAX, y_single, ORDER, PERIODIC, single_spls(iq))
        end do
        
        ! Test evaluation at multiple points
        do n_test = 1, 20
            x_eval = X_MIN + (X_MAX - X_MIN) * real(n_test-1, dp) / 19.0d0
            
            ! Batch evaluation
            call evaluate_batch_splines_1d(batch_spl, x_eval, y_batch_result)
            
            ! Compare with individual evaluations
            do iq = 1, N_QUANTITIES
                call evaluate_splines_1d(single_spls(iq), x_eval, y_single_result)
                
                if (abs(y_batch_result(iq) - y_single_result) > TOL) then
                    print *, "Evaluation mismatch at x=", x_eval, " iq=", iq
                    print *, "Batch: ", y_batch_result(iq)
                    print *, "Single: ", y_single_result
                    error stop "Evaluation mismatch"
                end if
            end do
        end do
        
        ! Clean up
        call destroy_batch_splines_1d(batch_spl)
        do iq = 1, N_QUANTITIES
            call destroy_splines_1d(single_spls(iq))
        end do
        
        print *, "  PASSED: 1D batch evaluation matches individual splines"
        
    end subroutine test_batch_spline_1d_evaluation
    
    
    subroutine test_batch_spline_1d_derivatives()
        use interpolate
        
        integer, parameter :: N_POINTS = 60
        integer, parameter :: N_QUANTITIES = 2
        integer, parameter :: ORDER = 5
        logical, parameter :: PERIODIC = .false.
        
        type(BatchSplineData1D) :: batch_spl
        type(SplineData1D) :: single_spls(N_QUANTITIES)
        
        real(dp) :: x(N_POINTS)
        real(dp) :: y_batch(N_POINTS, N_QUANTITIES)
        real(dp) :: y_single(N_POINTS)
        
        real(dp) :: x_eval
        real(dp) :: y_batch_result(N_QUANTITIES)
        real(dp) :: dy_batch_result(N_QUANTITIES)
        real(dp) :: d2y_batch_result(N_QUANTITIES)
        real(dp) :: y_single_result, dy_single_result, d2y_single_result
        
        integer :: iq, i, n_test
        
        print *, "Testing 1D batch spline derivatives..."
        
        ! Generate test data
        call linspace(X_MIN, X_MAX, N_POINTS, x)
        
        do iq = 1, N_QUANTITIES
            do i = 1, N_POINTS
                y_batch(i, iq) = cos(x(i)) * exp(-0.1d0*x(i)*real(iq, dp))
            end do
        end do
        
        ! Construct splines
        call construct_batch_splines_1d(X_MIN, X_MAX, y_batch, ORDER, PERIODIC, batch_spl)
        
        do iq = 1, N_QUANTITIES
            y_single(:) = y_batch(:, iq)
            call construct_splines_1d(X_MIN, X_MAX, y_single, ORDER, PERIODIC, single_spls(iq))
        end do
        
        ! Test derivatives at multiple points
        do n_test = 1, 15
            x_eval = X_MIN + (X_MAX - X_MIN) * real(n_test, dp) / 16.0d0
            
            ! Batch evaluation with derivatives
            call evaluate_batch_splines_1d_der(batch_spl, x_eval, y_batch_result, dy_batch_result)
            
            ! Compare first derivatives
            do iq = 1, N_QUANTITIES
                call evaluate_splines_1d_der(single_spls(iq), x_eval, y_single_result, dy_single_result)
                
                if (abs(y_batch_result(iq) - y_single_result) > TOL) then
                    error stop "Value mismatch in derivative test"
                end if
                if (abs(dy_batch_result(iq) - dy_single_result) > TOL) then
                    print *, "First derivative mismatch at x=", x_eval, " iq=", iq
                    print *, "Batch: ", dy_batch_result(iq)
                    print *, "Single: ", dy_single_result
                    error stop "First derivative mismatch"
                end if
            end do
            
            ! Test second derivatives
            call evaluate_batch_splines_1d_der2(batch_spl, x_eval, y_batch_result, dy_batch_result, d2y_batch_result)
            
            do iq = 1, N_QUANTITIES
                call evaluate_splines_1d_der2(single_spls(iq), x_eval, y_single_result, dy_single_result, d2y_single_result)
                
                if (abs(d2y_batch_result(iq) - d2y_single_result) > TOL) then
                    print *, "Second derivative mismatch at x=", x_eval, " iq=", iq
                    print *, "Batch: ", d2y_batch_result(iq)
                    print *, "Single: ", d2y_single_result
                    error stop "Second derivative mismatch"
                end if
            end do
        end do
        
        ! Clean up
        call destroy_batch_splines_1d(batch_spl)
        do iq = 1, N_QUANTITIES
            call destroy_splines_1d(single_spls(iq))
        end do
        
        print *, "  PASSED: 1D batch derivatives match individual splines"
        
    end subroutine test_batch_spline_1d_derivatives
    
    
    subroutine test_batch_spline_1d_single_extraction()
        use interpolate
        
        integer, parameter :: N_POINTS = 40
        integer, parameter :: N_QUANTITIES = 5
        integer, parameter :: ORDER = 3
        logical, parameter :: PERIODIC = .false.
        
        type(BatchSplineData1D) :: batch_spl
        type(SplineData1D) :: single_spls(N_QUANTITIES)
        
        real(dp) :: x(N_POINTS)
        real(dp) :: y_batch(N_POINTS, N_QUANTITIES)
        real(dp) :: y_single(N_POINTS)
        
        real(dp) :: x_eval
        real(dp) :: y_batch_all(N_QUANTITIES)
        real(dp) :: y_extracted, y_single_result
        
        integer :: iq, i, n_test, iq_test
        
        print *, "Testing 1D batch single quantity extraction..."
        
        ! Generate test data
        call linspace(X_MIN, X_MAX, N_POINTS, x)
        
        do iq = 1, N_QUANTITIES
            do i = 1, N_POINTS
                y_batch(i, iq) = x(i)**real(iq, dp) * exp(-x(i)/5.0d0)
            end do
        end do
        
        ! Construct splines
        call construct_batch_splines_1d(X_MIN, X_MAX, y_batch, ORDER, PERIODIC, batch_spl)
        
        do iq = 1, N_QUANTITIES
            y_single(:) = y_batch(:, iq)
            call construct_splines_1d(X_MIN, X_MAX, y_single, ORDER, PERIODIC, single_spls(iq))
        end do
        
        ! Test single extraction at multiple points
        do n_test = 1, 10
            x_eval = X_MIN + (X_MAX - X_MIN) * real(n_test, dp) / 11.0d0
            
            ! Full batch evaluation for reference
            call evaluate_batch_splines_1d(batch_spl, x_eval, y_batch_all)
            
            ! Test extraction of each quantity
            do iq_test = 1, N_QUANTITIES
                call evaluate_batch_splines_1d_single(batch_spl, x_eval, iq_test, y_extracted)
                
                ! Compare with full batch result
                if (abs(y_extracted - y_batch_all(iq_test)) > TOL) then
                    error stop "Single extraction doesn't match batch evaluation"
                end if
                
                ! Compare with individual spline
                call evaluate_splines_1d(single_spls(iq_test), x_eval, y_single_result)
                if (abs(y_extracted - y_single_result) > TOL) then
                    error stop "Single extraction doesn't match individual spline"
                end if
            end do
        end do
        
        ! Clean up
        call destroy_batch_splines_1d(batch_spl)
        do iq = 1, N_QUANTITIES
            call destroy_splines_1d(single_spls(iq))
        end do
        
        print *, "  PASSED: 1D batch single extraction works correctly"
        
    end subroutine test_batch_spline_1d_single_extraction
    
    
    subroutine test_batch_spline_1d_periodic()
        use interpolate
        
        integer, parameter :: N_POINTS = 45
        integer, parameter :: N_QUANTITIES = 3
        integer, parameter :: ORDER = 5
        logical, parameter :: PERIODIC = .true.
        
        type(BatchSplineData1D) :: batch_spl
        type(SplineData1D) :: single_spls(N_QUANTITIES)
        
        real(dp) :: x(N_POINTS)
        real(dp) :: y_batch(N_POINTS, N_QUANTITIES)
        real(dp) :: y_single(N_POINTS)
        
        real(dp) :: x_eval, x_eval_wrapped
        real(dp) :: y_batch_result(N_QUANTITIES)
        real(dp) :: y_batch_wrapped(N_QUANTITIES)
        
        integer :: iq, i
        real(dp) :: period
        
        print *, "Testing 1D batch spline with periodic boundary..."
        
        ! Generate periodic test data
        call linspace(0.0d0, TWOPI, N_POINTS, x)
        period = TWOPI
        
        do iq = 1, N_QUANTITIES
            do i = 1, N_POINTS
                y_batch(i, iq) = sin(x(i)*real(iq, dp))
            end do
        end do
        
        ! Make last point match first for periodicity
        do iq = 1, N_QUANTITIES
            y_batch(N_POINTS, iq) = y_batch(1, iq)
        end do
        
        ! Construct batch spline
        call construct_batch_splines_1d(0.0d0, TWOPI, y_batch, ORDER, PERIODIC, batch_spl)
        
        ! Test periodicity by evaluating outside domain
        x_eval = -0.5d0
        x_eval_wrapped = x_eval + period
        
        call evaluate_batch_splines_1d(batch_spl, x_eval, y_batch_result)
        call evaluate_batch_splines_1d(batch_spl, x_eval_wrapped, y_batch_wrapped)
        
        do iq = 1, N_QUANTITIES
            if (abs(y_batch_result(iq) - y_batch_wrapped(iq)) > TOL) then
                print *, "Periodic wrapping failed for iq=", iq
                print *, "At x=", x_eval, ": ", y_batch_result(iq)
                print *, "At wrapped x=", x_eval_wrapped, ": ", y_batch_wrapped(iq)
                error stop "Periodic boundary condition failed"
            end if
        end do
        
        ! Test at x = period + 0.3
        x_eval = period + 0.3d0
        x_eval_wrapped = 0.3d0
        
        call evaluate_batch_splines_1d(batch_spl, x_eval, y_batch_result)
        call evaluate_batch_splines_1d(batch_spl, x_eval_wrapped, y_batch_wrapped)
        
        do iq = 1, N_QUANTITIES
            if (abs(y_batch_result(iq) - y_batch_wrapped(iq)) > TOL) then
                error stop "Periodic boundary condition failed (upper)"
            end if
        end do
        
        ! Clean up
        call destroy_batch_splines_1d(batch_spl)
        
        print *, "  PASSED: 1D batch periodic boundary conditions work"
        
    end subroutine test_batch_spline_1d_periodic
    
    
    subroutine test_batch_spline_1d_memory_layout()
        use interpolate
        
        integer, parameter :: N_POINTS = 30
        integer, parameter :: N_QUANTITIES = 4
        integer, parameter :: ORDER = 3
        
        type(BatchSplineData1D) :: batch_spl
        real(dp) :: y_batch(N_POINTS, N_QUANTITIES)
        
        integer :: i, iq, k
        
        print *, "Testing 1D batch spline memory layout..."
        
        ! Generate test data
        do iq = 1, N_QUANTITIES
            do i = 1, N_POINTS
                y_batch(i, iq) = real(i*100 + iq, dp)
            end do
        end do
        
        ! Construct batch spline
        call construct_batch_splines_1d(X_MIN, X_MAX, y_batch, ORDER, .false., batch_spl)
        
        ! Verify allocation shape
        if (.not. allocated(batch_spl%coeff)) error stop "Coefficients not allocated"
        
        if (size(batch_spl%coeff, 1) /= N_QUANTITIES) error stop "Wrong coefficient quantities dimension"
        if (size(batch_spl%coeff, 2) /= ORDER + 1) error stop "Wrong coefficient order dimension"
        if (size(batch_spl%coeff, 3) /= N_POINTS) error stop "Wrong coefficient points dimension"
        
        ! The memory layout with quantities as the last dimension is automatically
        ! cache-friendly in Fortran due to column-major ordering.
        ! Accessing coefficients for all quantities at a given (k,i) point will be contiguous.
        
        ! Clean up
        call destroy_batch_splines_1d(batch_spl)
        
        print *, "  PASSED: 1D batch memory layout is cache-friendly"
        
    end subroutine test_batch_spline_1d_memory_layout
    
    
    subroutine bench_batch_vs_individual_1d()
        use interpolate
        
        integer, parameter :: N_POINTS = 100
        integer, parameter :: N_QUANTITIES = 6
        integer, parameter :: ORDER = 5
        integer, parameter :: N_EVALS = 10000
        
        type(BatchSplineData1D) :: batch_spl
        type(SplineData1D) :: single_spls(N_QUANTITIES)
        
        real(dp) :: x(N_POINTS)
        real(dp) :: y_batch(N_POINTS, N_QUANTITIES)
        real(dp) :: y_single(N_POINTS)
        
        real(dp) :: x_eval
        real(dp) :: y_batch_result(N_QUANTITIES)
        real(dp) :: y_single_result
        
        real(dp) :: t_start, t_end, t_batch, t_individual
        integer :: iq, i, n_eval
        real(dp) :: speedup
        
        print *, "Benchmarking batch vs individual 1D splines..."
        
        ! Generate test data
        call linspace(X_MIN, X_MAX, N_POINTS, x)
        
        do iq = 1, N_QUANTITIES
            do i = 1, N_POINTS
                y_batch(i, iq) = sin(x(i)*real(iq, dp)/2.0d0)
            end do
        end do
        
        ! Construct splines
        call construct_batch_splines_1d(X_MIN, X_MAX, y_batch, ORDER, .false., batch_spl)
        
        do iq = 1, N_QUANTITIES
            y_single(:) = y_batch(:, iq)
            call construct_splines_1d(X_MIN, X_MAX, y_single, ORDER, .false., single_spls(iq))
        end do
        
        ! Benchmark batch evaluation
        call cpu_time(t_start)
        do n_eval = 1, N_EVALS
            x_eval = X_MIN + mod(real(n_eval, dp)*0.01d0, X_MAX - X_MIN)
            call evaluate_batch_splines_1d(batch_spl, x_eval, y_batch_result)
        end do
        call cpu_time(t_end)
        t_batch = t_end - t_start
        
        ! Benchmark individual evaluations
        call cpu_time(t_start)
        do n_eval = 1, N_EVALS
            x_eval = X_MIN + mod(real(n_eval, dp)*0.01d0, X_MAX - X_MIN)
            do iq = 1, N_QUANTITIES
                call evaluate_splines_1d(single_spls(iq), x_eval, y_single_result)
            end do
        end do
        call cpu_time(t_end)
        t_individual = t_end - t_start
        
        speedup = t_individual / t_batch
        
        print *, "  Batch time: ", t_batch, " seconds"
        print *, "  Individual time: ", t_individual, " seconds"
        print *, "  Speedup: ", speedup, "x"
        
        if (speedup < 1.5d0) then
            print *, "  WARNING: Expected at least 1.5x speedup for 6 quantities"
        end if
        
        ! Clean up
        call destroy_batch_splines_1d(batch_spl)
        do iq = 1, N_QUANTITIES
            call destroy_splines_1d(single_spls(iq))
        end do
        
    end subroutine bench_batch_vs_individual_1d
    
    
    subroutine test_batch_spline_2d_construction()
        use interpolate
        
        integer, parameter :: N_POINTS(2) = [40, 45]
        integer, parameter :: N_QUANTITIES = 3
        integer, parameter :: ORDER(2) = [5, 3]
        logical, parameter :: PERIODIC(2) = [.false., .false.]
        
        type(BatchSplineData2D) :: batch_spl
        type(SplineData2D) :: single_spls(N_QUANTITIES)
        
        real(dp) :: y_batch(N_POINTS(1), N_POINTS(2), N_QUANTITIES)
        real(dp) :: y_single(N_POINTS(1), N_POINTS(2))
        
        real(dp) :: x1(N_POINTS(1)), x2(N_POINTS(2))
        integer :: iq, i1, i2, k1, k2
        
        print *, "Testing 2D batch spline construction..."
        
        ! Generate test data
        call linspace(X_MIN, X_MAX, N_POINTS(1), x1)
        call linspace(X_MIN, X_MAX, N_POINTS(2), x2)
        
        do iq = 1, N_QUANTITIES
            do i2 = 1, N_POINTS(2)
                do i1 = 1, N_POINTS(1)
                    y_batch(i1, i2, iq) = cos(x1(i1)*real(iq, dp)) * sin(x2(i2))
                end do
            end do
        end do
        
        ! Construct batch spline
        call construct_batch_splines_2d([X_MIN, X_MIN], [X_MAX, X_MAX], y_batch, ORDER, PERIODIC, batch_spl)
        
        ! Construct individual splines for comparison
        do iq = 1, N_QUANTITIES
            y_single(:,:) = y_batch(:,:,iq)
            call construct_splines_2d([X_MIN, X_MIN], [X_MAX, X_MAX], y_single, ORDER, PERIODIC, single_spls(iq))
        end do
        
        ! Verify metadata
        if (batch_spl%num_quantities /= N_QUANTITIES) error stop "Wrong num_quantities in 2D"
        if (any(batch_spl%num_points /= N_POINTS)) error stop "Wrong num_points in 2D"
        if (any(batch_spl%order /= ORDER)) error stop "Wrong order in 2D"
        
        ! Verify coefficients match
        do iq = 1, N_QUANTITIES
            do i2 = 1, N_POINTS(2)
                do i1 = 1, N_POINTS(1)
                    do k2 = 0, ORDER(2)
                        do k1 = 0, ORDER(1)
                            if (abs(batch_spl%coeff(iq, k1, k2, i1, i2) - &
                                   single_spls(iq)%coeff(k1, k2, i1, i2)) > TOL) then
                                error stop "2D coefficient mismatch"
                            end if
                        end do
                    end do
                end do
            end do
        end do
        
        ! Clean up
        call destroy_batch_splines_2d(batch_spl)
        do iq = 1, N_QUANTITIES
            call destroy_splines_2d(single_spls(iq))
        end do
        
        print *, "  PASSED: 2D batch construction matches individual splines"
        
    end subroutine test_batch_spline_2d_construction
    
    
    subroutine test_batch_spline_2d_evaluation()
        use interpolate
        
        integer, parameter :: N_POINTS(2) = [35, 40]
        integer, parameter :: N_QUANTITIES = 4
        integer, parameter :: ORDER(2) = [3, 3]
        
        type(BatchSplineData2D) :: batch_spl
        type(SplineData2D) :: single_spls(N_QUANTITIES)
        
        real(dp) :: y_batch(N_POINTS(1), N_POINTS(2), N_QUANTITIES)
        real(dp) :: y_single(N_POINTS(1), N_POINTS(2))
        
        real(dp) :: x1(N_POINTS(1)), x2(N_POINTS(2))
        real(dp) :: x_eval(2)
        real(dp) :: y_batch_result(N_QUANTITIES)
        real(dp) :: y_single_result
        
        integer :: iq, i1, i2, n_test
        
        print *, "Testing 2D batch spline evaluation..."
        
        ! Generate test data
        call linspace(X_MIN, X_MAX, N_POINTS(1), x1)
        call linspace(X_MIN, X_MAX, N_POINTS(2), x2)
        
        do iq = 1, N_QUANTITIES
            do i2 = 1, N_POINTS(2)
                do i1 = 1, N_POINTS(1)
                    y_batch(i1, i2, iq) = exp(-0.1d0*(x1(i1)+x2(i2))) * real(iq, dp)
                end do
            end do
        end do
        
        ! Construct splines
        call construct_batch_splines_2d([X_MIN, X_MIN], [X_MAX, X_MAX], y_batch, &
                                       ORDER, [.false., .false.], batch_spl)
        
        do iq = 1, N_QUANTITIES
            y_single(:,:) = y_batch(:,:,iq)
            call construct_splines_2d([X_MIN, X_MIN], [X_MAX, X_MAX], y_single, &
                                    ORDER, [.false., .false.], single_spls(iq))
        end do
        
        ! Test evaluation at multiple points
        do n_test = 1, 15
            x_eval(1) = X_MIN + (X_MAX - X_MIN) * real(n_test, dp) / 20.0d0
            x_eval(2) = X_MIN + (X_MAX - X_MIN) * real(n_test+3, dp) / 20.0d0
            
            ! Batch evaluation
            call evaluate_batch_splines_2d(batch_spl, x_eval, y_batch_result)
            
            ! Compare with individual evaluations
            do iq = 1, N_QUANTITIES
                call evaluate_splines_2d(single_spls(iq), x_eval, y_single_result)
                
                if (abs(y_batch_result(iq) - y_single_result) > TOL) then
                    error stop "2D evaluation mismatch"
                end if
            end do
        end do
        
        ! Clean up
        call destroy_batch_splines_2d(batch_spl)
        do iq = 1, N_QUANTITIES
            call destroy_splines_2d(single_spls(iq))
        end do
        
        print *, "  PASSED: 2D batch evaluation matches individual splines"
        
    end subroutine test_batch_spline_2d_evaluation
    
    
    subroutine test_batch_spline_2d_derivatives()
        use interpolate
        
        integer, parameter :: N_POINTS(2) = [30, 35]
        integer, parameter :: N_QUANTITIES = 2
        integer, parameter :: ORDER(2) = [5, 5]
        
        type(BatchSplineData2D) :: batch_spl
        type(SplineData2D) :: single_spls(N_QUANTITIES)
        
        real(dp) :: y_batch(N_POINTS(1), N_POINTS(2), N_QUANTITIES)
        real(dp) :: y_single(N_POINTS(1), N_POINTS(2))
        
        real(dp) :: x1(N_POINTS(1)), x2(N_POINTS(2))
        real(dp) :: x_eval(2)
        real(dp) :: y_batch_result(N_QUANTITIES)
        real(dp) :: dy_batch_result(2, N_QUANTITIES)
        real(dp) :: y_single_result
        real(dp) :: dy_single_result(2)
        
        integer :: iq, i1, i2, n_test
        
        print *, "Testing 2D batch spline derivatives..."
        
        ! Generate test data
        call linspace(X_MIN, X_MAX, N_POINTS(1), x1)
        call linspace(X_MIN, X_MAX, N_POINTS(2), x2)
        
        do iq = 1, N_QUANTITIES
            do i2 = 1, N_POINTS(2)
                do i1 = 1, N_POINTS(1)
                    y_batch(i1, i2, iq) = cos(x1(i1)) * sin(x2(i2)*real(iq, dp))
                end do
            end do
        end do
        
        ! Construct splines
        call construct_batch_splines_2d([X_MIN, X_MIN], [X_MAX, X_MAX], y_batch, &
                                       ORDER, [.false., .false.], batch_spl)
        
        do iq = 1, N_QUANTITIES
            y_single(:,:) = y_batch(:,:,iq)
            call construct_splines_2d([X_MIN, X_MIN], [X_MAX, X_MAX], y_single, &
                                    ORDER, [.false., .false.], single_spls(iq))
        end do
        
        ! Test derivatives at multiple points
        do n_test = 1, 10
            x_eval(1) = X_MIN + (X_MAX - X_MIN) * real(n_test+2, dp) / 15.0d0
            x_eval(2) = X_MIN + (X_MAX - X_MIN) * real(n_test+1, dp) / 15.0d0
            
            ! Batch evaluation with derivatives
            call evaluate_batch_splines_2d_der(batch_spl, x_eval, y_batch_result, dy_batch_result)
            
            ! Compare with individual evaluations
            do iq = 1, N_QUANTITIES
                call evaluate_splines_2d_der(single_spls(iq), x_eval, y_single_result, dy_single_result)
                
                if (abs(y_batch_result(iq) - y_single_result) > TOL) then
                    error stop "2D value mismatch in derivative test"
                end if
                
                if (any(abs(dy_batch_result(:, iq) - dy_single_result(:)) > TOL)) then
                    error stop "2D derivative mismatch"
                end if
            end do
        end do
        
        ! Clean up
        call destroy_batch_splines_2d(batch_spl)
        do iq = 1, N_QUANTITIES
            call destroy_splines_2d(single_spls(iq))
        end do
        
        print *, "  PASSED: 2D batch derivatives match individual splines"
        
    end subroutine test_batch_spline_2d_derivatives
    
    
    subroutine test_batch_spline_2d_mixed_periodic()
        use interpolate
        
        integer, parameter :: N_POINTS(2) = [30, 35]
        integer, parameter :: N_QUANTITIES = 2
        integer, parameter :: ORDER(2) = [3, 5]
        logical, parameter :: PERIODIC(2) = [.true., .false.]
        
        type(BatchSplineData2D) :: batch_spl
        
        real(dp) :: y_batch(N_POINTS(1), N_POINTS(2), N_QUANTITIES)
        real(dp) :: x1(N_POINTS(1)), x2(N_POINTS(2))
        real(dp) :: x_eval(2), x_eval_wrapped(2)
        real(dp) :: y_result(N_QUANTITIES), y_wrapped(N_QUANTITIES)
        
        integer :: iq, i1, i2
        
        print *, "Testing 2D batch spline with mixed periodic boundaries..."
        
        ! Generate test data (periodic in first dimension)
        call linspace(0.0d0, TWOPI, N_POINTS(1), x1)
        call linspace(X_MIN, X_MAX, N_POINTS(2), x2)
        
        do iq = 1, N_QUANTITIES
            do i2 = 1, N_POINTS(2)
                do i1 = 1, N_POINTS(1)
                    y_batch(i1, i2, iq) = sin(x1(i1)*real(iq, dp)) * exp(-x2(i2)/10.0d0)
                end do
            end do
            ! Enforce periodicity in first dimension
            y_batch(N_POINTS(1), :, iq) = y_batch(1, :, iq)
        end do
        
        ! Construct batch spline
        call construct_batch_splines_2d([0.0d0, X_MIN], [TWOPI, X_MAX], y_batch, ORDER, PERIODIC, batch_spl)
        
        ! Test periodic wrapping in first dimension
        x_eval = [TWOPI + 0.5d0, X_MIN + 1.0d0]
        x_eval_wrapped = [0.5d0, X_MIN + 1.0d0]
        
        call evaluate_batch_splines_2d(batch_spl, x_eval, y_result)
        call evaluate_batch_splines_2d(batch_spl, x_eval_wrapped, y_wrapped)
        
        do iq = 1, N_QUANTITIES
            if (abs(y_result(iq) - y_wrapped(iq)) > TOL) then
                error stop "Mixed periodic boundary failed"
            end if
        end do
        
        ! Clean up
        call destroy_batch_splines_2d(batch_spl)
        
        print *, "  PASSED: 2D batch mixed periodic boundaries work"
        
    end subroutine test_batch_spline_2d_mixed_periodic
    
    
    subroutine test_batch_spline_3d_construction()
        use interpolate
        
        integer, parameter :: N_POINTS(3) = [20, 25, 30]
        integer, parameter :: N_QUANTITIES = 3
        integer, parameter :: ORDER(3) = [3, 5, 3]
        logical, parameter :: PERIODIC(3) = [.false., .false., .false.]
        
        type(BatchSplineData3D) :: batch_spl
        type(SplineData3D) :: single_spls(N_QUANTITIES)
        
        real(dp) :: y_batch(N_POINTS(1), N_POINTS(2), N_POINTS(3), N_QUANTITIES)
        real(dp) :: y_single(N_POINTS(1), N_POINTS(2), N_POINTS(3))
        
        real(dp) :: x1(N_POINTS(1)), x2(N_POINTS(2)), x3(N_POINTS(3))
        integer :: iq, i1, i2, i3, k1, k2, k3
        
        print *, "Testing 3D batch spline construction..."
        
        ! Generate test data
        call linspace(X_MIN, X_MAX, N_POINTS(1), x1)
        call linspace(X_MIN, X_MAX, N_POINTS(2), x2)
        call linspace(X_MIN, X_MAX, N_POINTS(3), x3)
        
        do iq = 1, N_QUANTITIES
            do i3 = 1, N_POINTS(3)
                do i2 = 1, N_POINTS(2)
                    do i1 = 1, N_POINTS(1)
                        y_batch(i1, i2, i3, iq) = sin(x1(i1)) * cos(x2(i2)) * &
                                                  exp(-x3(i3)*real(iq, dp)/10.0d0)
                    end do
                end do
            end do
        end do
        
        ! Construct batch spline
        call construct_batch_splines_3d([X_MIN, X_MIN, X_MIN], [X_MAX, X_MAX, X_MAX], &
                                       y_batch, ORDER, PERIODIC, batch_spl)
        
        ! Verify metadata
        if (batch_spl%num_quantities /= N_QUANTITIES) error stop "Wrong num_quantities in 3D"
        if (any(batch_spl%num_points /= N_POINTS)) error stop "Wrong num_points in 3D"
        if (any(batch_spl%order /= ORDER)) error stop "Wrong order in 3D"
        
        ! Clean up
        call destroy_batch_splines_3d(batch_spl)
        
        print *, "  PASSED: 3D batch construction completed successfully"
        
    end subroutine test_batch_spline_3d_construction
    
    
    subroutine test_batch_spline_3d_evaluation()
        use interpolate
        
        integer, parameter :: N_POINTS(3) = [25, 20, 22]
        integer, parameter :: N_QUANTITIES = 2
        integer, parameter :: ORDER(3) = [3, 3, 3]
        
        type(BatchSplineData3D) :: batch_spl
        type(SplineData3D) :: single_spls(N_QUANTITIES)
        
        real(dp) :: y_batch(N_POINTS(1), N_POINTS(2), N_POINTS(3), N_QUANTITIES)
        real(dp) :: y_single(N_POINTS(1), N_POINTS(2), N_POINTS(3))
        
        real(dp) :: x1(N_POINTS(1)), x2(N_POINTS(2)), x3(N_POINTS(3))
        real(dp) :: x_eval(3)
        real(dp) :: y_batch_result(N_QUANTITIES)
        real(dp) :: y_single_result
        
        integer :: iq, i1, i2, i3, n_test
        
        print *, "Testing 3D batch spline evaluation..."
        
        ! Generate test data
        call linspace(X_MIN, X_MAX, N_POINTS(1), x1)
        call linspace(X_MIN, X_MAX, N_POINTS(2), x2)
        call linspace(X_MIN, X_MAX, N_POINTS(3), x3)
        
        do iq = 1, N_QUANTITIES
            do i3 = 1, N_POINTS(3)
                do i2 = 1, N_POINTS(2)
                    do i1 = 1, N_POINTS(1)
                        y_batch(i1, i2, i3, iq) = x1(i1)**real(iq, dp) + x2(i2) + x3(i3)
                    end do
                end do
            end do
        end do
        
        ! Construct splines
        call construct_batch_splines_3d([X_MIN, X_MIN, X_MIN], [X_MAX, X_MAX, X_MAX], &
                                       y_batch, ORDER, [.false., .false., .false.], batch_spl)
        
        do iq = 1, N_QUANTITIES
            y_single(:,:,:) = y_batch(:,:,:,iq)
            call construct_splines_3d([X_MIN, X_MIN, X_MIN], [X_MAX, X_MAX, X_MAX], &
                                    y_single, ORDER, [.false., .false., .false.], single_spls(iq))
        end do
        
        ! Test evaluation at multiple points
        do n_test = 1, 10
            x_eval(1) = X_MIN + (X_MAX - X_MIN) * real(n_test, dp) / 12.0d0
            x_eval(2) = X_MIN + (X_MAX - X_MIN) * real(n_test+1, dp) / 12.0d0
            x_eval(3) = X_MIN + (X_MAX - X_MIN) * real(n_test+2, dp) / 12.0d0
            
            ! Batch evaluation
            call evaluate_batch_splines_3d(batch_spl, x_eval, y_batch_result)
            
            ! Compare with individual evaluations
            do iq = 1, N_QUANTITIES
                call evaluate_splines_3d(single_spls(iq), x_eval, y_single_result)
                
                if (abs(y_batch_result(iq) - y_single_result) > TOL) then
                    error stop "3D evaluation mismatch"
                end if
            end do
        end do
        
        ! Clean up
        call destroy_batch_splines_3d(batch_spl)
        do iq = 1, N_QUANTITIES
            call destroy_splines_3d(single_spls(iq))
        end do
        
        print *, "  PASSED: 3D batch evaluation matches individual splines"
        
    end subroutine test_batch_spline_3d_evaluation
    
    
    subroutine test_batch_spline_3d_derivatives()
        use interpolate
        
        integer, parameter :: N_POINTS(3) = [20, 22, 24]
        integer, parameter :: N_QUANTITIES = 2
        integer, parameter :: ORDER(3) = [5, 3, 5]
        
        type(BatchSplineData3D) :: batch_spl
        type(SplineData3D) :: single_spls(N_QUANTITIES)
        
        real(dp) :: y_batch(N_POINTS(1), N_POINTS(2), N_POINTS(3), N_QUANTITIES)
        real(dp) :: y_single(N_POINTS(1), N_POINTS(2), N_POINTS(3))
        
        real(dp) :: x1(N_POINTS(1)), x2(N_POINTS(2)), x3(N_POINTS(3))
        real(dp) :: x_eval(3)
        real(dp) :: y_batch_result(N_QUANTITIES)
        real(dp) :: dy_batch_result(3, N_QUANTITIES)
        real(dp) :: d2y_batch_result(6, N_QUANTITIES)
        real(dp) :: y_single_result
        real(dp) :: dy_single_result(3)
        real(dp) :: d2y_single_result(6)
        
        integer :: iq, i1, i2, i3, n_test
        
        print *, "Testing 3D batch spline derivatives..."
        
        ! Generate test data
        call linspace(X_MIN, X_MAX, N_POINTS(1), x1)
        call linspace(X_MIN, X_MAX, N_POINTS(2), x2)
        call linspace(X_MIN, X_MAX, N_POINTS(3), x3)
        
        do iq = 1, N_QUANTITIES
            do i3 = 1, N_POINTS(3)
                do i2 = 1, N_POINTS(2)
                    do i1 = 1, N_POINTS(1)
                        y_batch(i1, i2, i3, iq) = cos(x1(i1)*real(iq, dp)) * &
                                                  sin(x2(i2)) * cos(x3(i3))
                    end do
                end do
            end do
        end do
        
        ! Construct splines
        call construct_batch_splines_3d([X_MIN, X_MIN, X_MIN], [X_MAX, X_MAX, X_MAX], &
                                       y_batch, ORDER, [.false., .false., .false.], batch_spl)
        
        do iq = 1, N_QUANTITIES
            y_single(:,:,:) = y_batch(:,:,:,iq)
            call construct_splines_3d([X_MIN, X_MIN, X_MIN], [X_MAX, X_MAX, X_MAX], &
                                    y_single, ORDER, [.false., .false., .false.], single_spls(iq))
        end do
        
        ! Test derivatives
        do n_test = 1, 5
            x_eval(1) = X_MIN + (X_MAX - X_MIN) * real(n_test+1, dp) / 8.0d0
            x_eval(2) = X_MIN + (X_MAX - X_MIN) * real(n_test+2, dp) / 8.0d0
            x_eval(3) = X_MIN + (X_MAX - X_MIN) * real(n_test, dp) / 8.0d0
            
            ! Batch evaluation with derivatives
            call evaluate_batch_splines_3d_der(batch_spl, x_eval, y_batch_result, &
                                               dy_batch_result)
            
            ! Compare with individual evaluations
            do iq = 1, N_QUANTITIES
                call evaluate_splines_3d_der(single_spls(iq), x_eval, y_single_result, &
                                             dy_single_result)
                
                if (abs(y_batch_result(iq) - y_single_result) > TOL) then
                    error stop "3D value mismatch in derivative test"
                end if
                
                if (any(abs(dy_batch_result(:, iq) - dy_single_result(:)) > TOL)) then
                    error stop "3D first derivative mismatch"
                end if
            end do
        end do
        
        ! Clean up
        call destroy_batch_splines_3d(batch_spl)
        do iq = 1, N_QUANTITIES
            call destroy_splines_3d(single_spls(iq))
        end do
        
        print *, "  PASSED: 3D batch derivatives match individual splines"
        
    end subroutine test_batch_spline_3d_derivatives
    
    
    subroutine test_batch_spline_3d_field_components()
        use interpolate
        
        ! Simulate magnetic field with A and B components
        integer, parameter :: N_POINTS(3) = [30, 32, 35]
        integer, parameter :: ORDER(3) = [5, 5, 5]
        logical, parameter :: PERIODIC(3) = [.false., .true., .false.]
        
        type(BatchSplineData3D) :: A_batch_spl, B_batch_spl
        
        real(dp) :: A_batch(N_POINTS(1), N_POINTS(2), N_POINTS(3), 3)
        real(dp) :: B_batch(N_POINTS(1), N_POINTS(2), N_POINTS(3), 3)
        
        real(dp) :: x1(N_POINTS(1)), x2(N_POINTS(2)), x3(N_POINTS(3))
        real(dp) :: x_eval(3)
        real(dp) :: A_result(3), B_result(3)
        
        integer :: i1, i2, i3
        real(dp) :: r, theta, z
        
        print *, "Testing 3D batch splines for field components..."
        
        ! Generate cylindrical-like grid
        call linspace(1.0d0, 2.0d0, N_POINTS(1), x1)  ! r
        call linspace(0.0d0, TWOPI, N_POINTS(2), x2)  ! theta (periodic)
        call linspace(-1.0d0, 1.0d0, N_POINTS(3), x3) ! z
        
        ! Generate field components
        do i3 = 1, N_POINTS(3)
            z = x3(i3)
            do i2 = 1, N_POINTS(2)
                theta = x2(i2)
                do i1 = 1, N_POINTS(1)
                    r = x1(i1)
                    
                    ! Vector potential A
                    A_batch(i1, i2, i3, 1) = 0.0d0  ! A_r
                    A_batch(i1, i2, i3, 2) = -z/r   ! A_theta
                    A_batch(i1, i2, i3, 3) = log(r) ! A_z
                    
                    ! Magnetic field B
                    B_batch(i1, i2, i3, 1) = 1.0d0/r        ! B_r
                    B_batch(i1, i2, i3, 2) = 0.0d0          ! B_theta
                    B_batch(i1, i2, i3, 3) = 1.0d0/(r**2)   ! B_z
                end do
            end do
        end do
        
        ! Enforce periodicity in theta
        A_batch(:, N_POINTS(2), :, :) = A_batch(:, 1, :, :)
        B_batch(:, N_POINTS(2), :, :) = B_batch(:, 1, :, :)
        
        ! Construct batch splines
        call construct_batch_splines_3d([1.0d0, 0.0d0, -1.0d0], [2.0d0, TWOPI, 1.0d0], &
                                       A_batch, ORDER, PERIODIC, A_batch_spl)
        call construct_batch_splines_3d([1.0d0, 0.0d0, -1.0d0], [2.0d0, TWOPI, 1.0d0], &
                                       B_batch, ORDER, PERIODIC, B_batch_spl)
        
        ! Evaluate at test point
        x_eval = [1.5d0, PI, 0.3d0]
        
        call evaluate_batch_splines_3d(A_batch_spl, x_eval, A_result)
        call evaluate_batch_splines_3d(B_batch_spl, x_eval, B_result)
        
        ! Verify reasonable values
        if (any(abs(A_result) > 10.0d0)) error stop "Unreasonable A field values"
        if (any(abs(B_result) > 10.0d0)) error stop "Unreasonable B field values"
        
        ! Test periodic wrapping in theta
        x_eval(2) = TWOPI + 0.1d0  ! Should wrap to 0.1
        call evaluate_batch_splines_3d(A_batch_spl, x_eval, A_result)
        
        x_eval(2) = 0.1d0
        call evaluate_batch_splines_3d(A_batch_spl, x_eval, B_result)
        
        if (any(abs(A_result - B_result) > TOL)) then
            error stop "Periodic wrapping failed for field components"
        end if
        
        ! Clean up
        call destroy_batch_splines_3d(A_batch_spl)
        call destroy_batch_splines_3d(B_batch_spl)
        
        print *, "  PASSED: 3D batch field components work correctly"
        
    end subroutine test_batch_spline_3d_field_components

end program test_batch_interpolate