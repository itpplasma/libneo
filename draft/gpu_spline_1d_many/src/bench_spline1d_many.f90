program bench_spline1d_many
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer
    use batch_interpolate_types, only: BatchSplineData1D
    use batch_interpolate_1d, only: construct_batch_splines_1d, &
                                    construct_batch_splines_1d_lines, &
                                    construct_batch_splines_1d_resident, &
                                    construct_batch_splines_1d_resident_device, &
                                    destroy_batch_splines_1d
    use draft_batch_splines_many_api, only: evaluate_batch_splines_1d_many
    use spline1d_many_offload, only: spline1d_many_setup, &
                                     spline1d_many_teardown, &
                                     spline1d_many_eval_host, &
                                     spline1d_many_eval_resident
    use util_lcg_rng, only: lcg_state_t, lcg_init, lcg_uniform_0_1
    use util_timer, only: wall_time
    use util_bench_args, only: get_env_int
    implicit none

    integer, parameter :: num_points = 2048
    integer, parameter :: npts = 2000000
    integer, parameter :: niter = 20
    integer, parameter :: nbuild = 3
    integer, parameter :: nbuild_resident = 3
    integer, parameter :: nbuild_repeat = 5
    logical, parameter :: periodic = .true.

    real(dp), parameter :: x_min = 1.23d0
    real(dp), parameter :: x_max = 40.0d0 + x_min

    type(BatchSplineData1D) :: spl
    integer :: order
    real(dp), allocatable :: y_batch(:)
    real(dp), allocatable :: x_eval(:)
    real(dp), allocatable, target :: y_ref(:)
    real(dp), pointer :: y_ref2d(:, :)
    type(c_ptr) :: y_ref_ptr

    real(dp), allocatable :: x_grid(:)
    real(dp), allocatable :: y_grid(:, :)

    type(lcg_state_t) :: rng
    integer :: ip, iq
    integer :: num_quantities
    real(dp) :: best
    real(dp) :: diff_max

    order = 5
    call parse_order_arg(order)
    num_quantities = get_env_int("LIBNEO_BENCH_NQ", 8)
    if (num_quantities < 1) then
        error stop "bench_spline1d_many: LIBNEO_BENCH_NQ must be >= 1"
    end if

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

    call bench_build_1d(x_min, x_max, y_grid, spl)
    call bench_build_lines_1d(x_min, x_max, y_grid)
    call bench_build_resident_1d(x_min, x_max, y_grid, spl)
    call bench_build_device_1d()

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
#endif

    call destroy_batch_splines_1d(spl)

contains

    subroutine parse_order_arg(order_out)
        integer, intent(inout) :: order_out

        character(len=32) :: arg
        integer :: stat

        call get_command_argument(1, arg)
        if (len_trim(arg) == 0) return

        read (arg, *, iostat=stat) order_out
        if (stat /= 0) then
            error stop "bench_spline1d_many: failed to parse order argument"
        end if
        if (order_out < 3 .or. order_out > 5) then
            error stop "bench_spline1d_many: order must be 3..5"
        end if
    end subroutine parse_order_arg

    subroutine bench_build_1d(x_min_local, x_max_local, y_grid_local, spl_out)
        real(dp), intent(in) :: x_min_local
        real(dp), intent(in) :: x_max_local
        real(dp), intent(in) :: y_grid_local(:, :)
        type(BatchSplineData1D), intent(out) :: spl_out

        integer :: it, rep
        real(dp) :: t0, t1, dt, best_build
        real(dp) :: grid_pts_per_s

        best_build = huge(1.0d0)
        do it = 1, nbuild
            t0 = wall_time()
            do rep = 1, nbuild_repeat
                call construct_batch_splines_1d(x_min_local, x_max_local, &
                                                y_grid_local, order, periodic, spl_out)
                call destroy_batch_splines_1d(spl_out)
            end do
            t1 = wall_time()
            dt = (t1 - t0) / real(nbuild_repeat, dp)
            best_build = min(best_build, dt)
        end do

        call construct_batch_splines_1d(x_min_local, x_max_local, y_grid_local, &
                                        order, periodic, spl_out)

        grid_pts_per_s = real(num_points*num_quantities, dp) / best_build
        print *, "build best_s ", best_build, " grid_pts_per_s ", grid_pts_per_s
    end subroutine bench_build_1d

    subroutine bench_build_lines_1d(x_min_local, x_max_local, y_grid_local)
        real(dp), intent(in) :: x_min_local
        real(dp), intent(in) :: x_max_local
        real(dp), intent(in) :: y_grid_local(:, :)

        type(BatchSplineData1D) :: spl_tmp
        integer :: it, rep
        real(dp) :: t0, t1, dt, best_build
        real(dp) :: grid_pts_per_s

        best_build = huge(1.0d0)
        do it = 1, nbuild
            t0 = wall_time()
            do rep = 1, nbuild_repeat
                call construct_batch_splines_1d_lines(x_min_local, x_max_local, &
                                                      y_grid_local, order, periodic, &
                                                      spl_tmp)
                call destroy_batch_splines_1d(spl_tmp)
            end do
            t1 = wall_time()
            dt = (t1 - t0) / real(nbuild_repeat, dp)
            best_build = min(best_build, dt)
        end do

        grid_pts_per_s = real(num_points*num_quantities, dp) / best_build
        print *, "build_lines best_s ", best_build, " grid_pts_per_s ", grid_pts_per_s
    end subroutine bench_build_lines_1d

    subroutine bench_build_resident_1d(x_min_local, x_max_local, y_grid_local, spl_out)
        real(dp), intent(in) :: x_min_local
        real(dp), intent(in) :: x_max_local
        real(dp), intent(in) :: y_grid_local(:, :)
        type(BatchSplineData1D), intent(inout) :: spl_out

        integer :: it, rep
        real(dp) :: t0, t1, dt, best_build
        real(dp) :: grid_pts_per_s

#if defined(LIBNEO_ENABLE_OPENACC)
        call destroy_batch_splines_1d(spl_out)

        best_build = huge(1.0d0)
        do it = 1, nbuild_resident
            t0 = wall_time()
            do rep = 1, nbuild_repeat
                call construct_batch_splines_1d_resident(x_min_local, x_max_local, &
                                                         y_grid_local, &
                                                         order, periodic, spl_out)
                !$acc wait
                call destroy_batch_splines_1d(spl_out)
            end do
            t1 = wall_time()
            dt = (t1 - t0) / real(nbuild_repeat, dp)
            best_build = min(best_build, dt)
        end do

        call construct_batch_splines_1d_resident(x_min_local, x_max_local, &
                                                 y_grid_local, order, periodic, spl_out)

        grid_pts_per_s = real(num_points*num_quantities, dp) / best_build
        print *, "build_resident best_s ", best_build, " grid_pts_per_s ", &
            grid_pts_per_s
#else
        call destroy_batch_splines_1d(spl_out)
        call construct_batch_splines_1d(x_min_local, x_max_local, y_grid_local, &
                                        order, periodic, spl_out)
#endif
    end subroutine bench_build_resident_1d

    subroutine bench_build_device_1d()
        type(BatchSplineData1D), save :: spl_dev
        integer :: it, rep
        real(dp) :: t0, t1, dt, best_build
        real(dp) :: grid_pts_per_s

#if defined(LIBNEO_ENABLE_OPENACC)
        best_build = huge(1.0d0)
        do it = 1, nbuild_resident
            t0 = wall_time()
            do rep = 1, nbuild_repeat
                call construct_batch_splines_1d_resident_device(&
                    x_min, x_max, y_grid, order, periodic, spl_dev, &
                    update_host=.true., assume_y_present=.false.)
                !$acc wait
                call destroy_batch_splines_1d(spl_dev)
            end do
            t1 = wall_time()
            dt = (t1 - t0) / real(nbuild_repeat, dp)
            best_build = min(best_build, dt)
        end do
        grid_pts_per_s = real(num_points*num_quantities, dp) / best_build
        print *, "build_device_host best_s ", best_build, " grid_pts_per_s ", &
            grid_pts_per_s

        !$acc enter data copyin(y_grid)
        best_build = huge(1.0d0)
        do it = 1, nbuild_resident
            t0 = wall_time()
            do rep = 1, nbuild_repeat
                call construct_batch_splines_1d_resident_device(&
                    x_min, x_max, y_grid, order, periodic, spl_dev, &
                    update_host=.false., assume_y_present=.true.)
                !$acc wait
                call destroy_batch_splines_1d(spl_dev)
            end do
            t1 = wall_time()
            dt = (t1 - t0) / real(nbuild_repeat, dp)
            best_build = min(best_build, dt)
        end do
        !$acc exit data delete(y_grid)
        grid_pts_per_s = real(num_points*num_quantities, dp) / best_build
        print *, "build_device_gpu best_s ", best_build, " grid_pts_per_s ", &
            grid_pts_per_s
#endif
    end subroutine bench_build_device_1d

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
                                         spl_local%x_min, spl_local%h_step, &
                                         spl_local%coeff, &
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
        real(dp) :: tsetup0, tsetup1

        tsetup0 = wall_time()
        call spline1d_many_setup(spl_local%coeff, x, y)
        !$acc wait
        tsetup1 = wall_time()
        print *, "openacc setup_s ", tsetup1 - tsetup0
        call spline1d_many_eval_resident(spl_local%order, spl_local%num_points, &
                                         spl_local%num_quantities, spl_local%periodic, &
                                         spl_local%x_min, spl_local%h_step, &
                                         spl_local%coeff, &
                                         x, y)
        !$acc wait

        best = huge(1.0d0)
        do it = 1, niter
            tstart = wall_time()
            call spline1d_many_eval_resident(spl_local%order, spl_local%num_points, &
                                             spl_local%num_quantities, &
                                             spl_local%periodic, &
                                             spl_local%x_min, spl_local%h_step, &
                                             spl_local%coeff, &
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

end program bench_spline1d_many
