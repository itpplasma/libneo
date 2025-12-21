program bench_spline1d_many
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer
    use batch_interpolate_types, only: BatchSplineData1D
    use batch_interpolate_1d, only: construct_batch_splines_1d, destroy_batch_splines_1d
    use draft_batch_splines_many_api, only: evaluate_batch_splines_1d_many
    use spline1d_many_offload, only: spline1d_many_setup, &
                                     spline1d_many_teardown, &
                                     spline1d_many_eval_host, &
                                     spline1d_many_eval_resident
    use util_lcg_rng, only: lcg_state_t, lcg_init, lcg_uniform_0_1
    implicit none

    integer, parameter :: order = 5
    integer, parameter :: num_points = 2048
    integer, parameter :: num_quantities = 8
    integer, parameter :: npts = 2000000
    integer, parameter :: niter = 20
    logical, parameter :: periodic = .true.

    real(dp), parameter :: x_min = 1.23d0
    real(dp), parameter :: x_max = 40.0d0 + x_min

    type(BatchSplineData1D) :: spl
    real(dp), allocatable :: y_batch(:)
    real(dp), allocatable :: x_eval(:)
    real(dp), allocatable, target :: y_ref(:)
    real(dp), pointer :: y_ref2d(:, :)
    type(c_ptr) :: y_ref_ptr

    real(dp), allocatable :: x_grid(:)
    real(dp), allocatable :: y_grid(:, :)

    type(lcg_state_t) :: rng
    integer :: ip, iq
    real(dp) :: best
    real(dp) :: diff_max

    allocate (x_grid(num_points))
    allocate (y_grid(num_points, num_quantities))

    do ip = 1, num_points
        x_grid(ip) = x_min + (x_max - x_min)*real(ip - 1, dp)/real(num_points - 1, dp)
    end do
    do iq = 1, num_quantities
        do ip = 1, num_points
            y_grid(ip, iq) = cos(x_grid(ip) + real(iq - 1, dp)*0.1d0)
        end do
    end do

    call construct_batch_splines_1d(x_min, x_max, y_grid, order, periodic, spl)

    allocate (x_eval(npts))
    allocate (y_ref(num_quantities*npts))
    allocate (y_batch(num_quantities*npts))

    call lcg_init(rng, 1234567_int64)
    do ip = 1, npts
        x_eval(ip) = x_min + (x_max - x_min)*lcg_uniform_0_1(rng)
    end do

    y_ref_ptr = c_loc(y_ref(1))
    call c_f_pointer(y_ref_ptr, y_ref2d, [spl%num_quantities, npts])
    call evaluate_batch_splines_1d_many(spl, x_eval, y_ref2d)

    print *, "Benchmark parameters"
    print *, "order         ", order
    print *, "num_points    ", num_points
    print *, "num_quantities", num_quantities
    print *, "npts          ", npts
    print *, "niter         ", niter
    print *, "periodic      ", periodic

    call bench_cpu(spl, x_eval, y_batch, y_ref)
#if defined(LIBNEO_ENABLE_OPENACC)
    call bench_openacc(spl, x_eval, y_batch, y_ref)
#elif defined(LIBNEO_ENABLE_OPENMP)
    call bench_openmp(spl, x_eval, y_batch, y_ref)
#endif

    call destroy_batch_splines_1d(spl)

contains

    real(dp) function wall_time() result(t)
        integer :: count, rate, max_count
        call system_clock(count, rate, max_count)
        t = real(count, dp)/real(rate, dp)
    end function wall_time

    subroutine bench_cpu(spl_local, x, y, y_expected)
        type(BatchSplineData1D), intent(in) :: spl_local
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(:)
        real(dp), intent(in) :: y_expected(:)

        integer :: it
        real(dp) :: tstart, tend, dt

        best = huge(1.0d0)
        do it = 1, niter
            tstart = wall_time()
            call spline1d_many_eval_host(spl_local%order, spl_local%num_points, &
                                         spl_local%num_quantities, spl_local%periodic, &
                                         spl_local%x_min, spl_local%h_step, spl_local%coeff, &
                                         x, y)
            tend = wall_time()
            dt = tend - tstart
            best = min(best, dt)
        end do

        diff_max = maxval(abs(y - y_expected))
        print *, "cpu best_s ", best, " pts_per_s ", real(size(x), dp)/best, &
            " max_abs_diff ", diff_max
    end subroutine bench_cpu

    subroutine bench_openacc(spl_local, x, y, y_expected)
        type(BatchSplineData1D), intent(in) :: spl_local
        real(dp), intent(in), target :: x(:)
        real(dp), intent(inout), target :: y(:)
        real(dp), intent(in) :: y_expected(:)

        integer :: it
        real(dp) :: tstart, tend, dt

        call spline1d_many_setup(spl_local%coeff, x, y)
        call spline1d_many_eval_resident(spl_local%order, spl_local%num_points, &
                                         spl_local%num_quantities, spl_local%periodic, &
                                         spl_local%x_min, spl_local%h_step, spl_local%coeff, &
                                         x, y)
        !$acc wait

        best = huge(1.0d0)
        do it = 1, niter
            tstart = wall_time()
            call spline1d_many_eval_resident(spl_local%order, spl_local%num_points, &
                                             spl_local%num_quantities, spl_local%periodic, &
                                             spl_local%x_min, spl_local%h_step, spl_local%coeff, &
                                             x, y)
            !$acc wait
            tend = wall_time()
            dt = tend - tstart
            best = min(best, dt)
        end do

        !$acc update self(y)
        diff_max = maxval(abs(y - y_expected))
        print *, "openacc best_s ", best, " pts_per_s ", real(size(x), dp)/best, &
            " max_abs_diff ", diff_max

        call spline1d_many_teardown(spl_local%coeff, x, y)
    end subroutine bench_openacc

    subroutine bench_openmp(spl_local, x, y, y_expected)
        type(BatchSplineData1D), intent(in) :: spl_local
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(:)
        real(dp), intent(in) :: y_expected(:)

        integer :: it
        real(dp) :: tstart, tend, dt

        call spline1d_many_setup(spl_local%coeff, x, y)
        call spline1d_many_eval_resident(spl_local%order, spl_local%num_points, &
                                         spl_local%num_quantities, spl_local%periodic, &
                                         spl_local%x_min, spl_local%h_step, spl_local%coeff, &
                                         x, y)

        best = huge(1.0d0)
        do it = 1, niter
            tstart = wall_time()
            call spline1d_many_eval_resident(spl_local%order, spl_local%num_points, &
                                             spl_local%num_quantities, spl_local%periodic, &
                                             spl_local%x_min, spl_local%h_step, spl_local%coeff, &
                                             x, y)
            tend = wall_time()
            dt = tend - tstart
            best = min(best, dt)
        end do

        !$omp target update from(y)
        diff_max = maxval(abs(y - y_expected))
        print *, "openmp best_s ", best, " pts_per_s ", real(size(x), dp)/best, &
            " max_abs_diff ", diff_max

        call spline1d_many_teardown(spl_local%coeff, x, y)
    end subroutine bench_openmp

end program bench_spline1d_many
