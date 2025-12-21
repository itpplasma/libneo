program bench_spline1d_many
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer
    use batch_interpolate_types, only: BatchSplineData1D
    use batch_interpolate_1d, only: construct_batch_splines_1d, &
                                    construct_batch_splines_1d_lines, &
                                    construct_batch_splines_1d_resident, &
                                    construct_batch_splines_1d_resident_device, &
                                    evaluate_batch_splines_1d, &
                                    evaluate_batch_splines_1d_many, &
                                    destroy_batch_splines_1d
    use spline1d_many_offload, only: spline1d_many_setup, &
                                     spline1d_many_teardown, &
                                     spline1d_many_eval_host, &
                                     spline1d_many_eval_resident
    use util_lcg_rng, only: lcg_state_t, lcg_init, lcg_uniform_0_1
    use util_timer, only: wall_time
    use util_bench_args, only: get_env_int
    implicit none

    integer :: num_points
    integer :: npts
    integer :: niter
    integer :: nbuild
    integer :: nbuild_resident
    integer :: nbuild_repeat
    logical, parameter :: periodic = .true.

    real(dp), parameter :: x_min = 1.23d0
    real(dp), parameter :: x_max = 40.0d0 + x_min

    type(BatchSplineData1D) :: spl
    integer :: order
    real(dp), allocatable, target :: y_batch(:)
    real(dp), allocatable, target :: y_public(:)
    real(dp), allocatable :: x_eval(:)
    real(dp), allocatable, target :: y_ref(:)
    type(c_ptr) :: y_ptr

    real(dp), allocatable :: x_grid(:)
    real(dp), allocatable :: y_grid(:, :)

    type(lcg_state_t) :: rng
    integer :: ip, iq
    integer :: num_quantities
    integer :: npts_old
    real(dp) :: best
    real(dp) :: diff_max

    num_points = get_env_int("LIBNEO_BENCH_NUM_POINTS", 2048)
    if (num_points < 8) error stop "bench_spline1d_many: LIBNEO_BENCH_NUM_POINTS < 8"

    npts = get_env_int("LIBNEO_BENCH_NPTS", 2000000)
    if (npts < 1) error stop "bench_spline1d_many: LIBNEO_BENCH_NPTS < 1"

    niter = get_env_int("LIBNEO_BENCH_NITER", 20)
    if (niter < 1) error stop "bench_spline1d_many: LIBNEO_BENCH_NITER < 1"

    nbuild = get_env_int("LIBNEO_BENCH_NBUILD", 3)
    if (nbuild < 1) error stop "bench_spline1d_many: LIBNEO_BENCH_NBUILD < 1"

    nbuild_resident = get_env_int("LIBNEO_BENCH_NBUILD_RESIDENT", nbuild)
    if (nbuild_resident < 1) error stop "bench_spline1d_many: LIBNEO_BENCH_NBUILD_RESIDENT < 1"

    nbuild_repeat = get_env_int("LIBNEO_BENCH_NBUILD_REPEAT", 5)
    if (nbuild_repeat < 1) error stop "bench_spline1d_many: LIBNEO_BENCH_NBUILD_REPEAT < 1"

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
    allocate (y_public(num_quantities*npts))

    call lcg_init(rng, 1234567_int64)
    do ip = 1, npts
        x_eval(ip) = x_min + (x_max - x_min)*lcg_uniform_0_1(rng)
    end do

    call spline1d_many_eval_host(spl%order, spl%num_points, spl%num_quantities, &
                                 spl%periodic, spl%x_min, spl%h_step, spl%coeff, &
                                 x_eval, y_ref)

    print *, "Benchmark parameters"
    print *, "order         ", order
    print *, "num_points    ", num_points
    print *, "num_quantities", num_quantities
    print *, "npts          ", npts
    print *, "niter         ", niter
    print *, "periodic      ", periodic

    npts_old = min(npts, get_env_int("LIBNEO_BENCH_OLD_NPTS", 200000))
    print *, "npts_old      ", npts_old

    call bench_cpu(spl, x_eval, y_batch, y_ref)
    call bench_old_public(spl, x_eval, npts_old)
#if defined(LIBNEO_ENABLE_OPENACC)
    call bench_openacc(spl, x_eval, y_batch, y_ref)
#endif
    call bench_public(spl, x_eval, y_public, y_ref)

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

    subroutine bench_old_public(spl_local, x, npts_local)
        type(BatchSplineData1D), intent(in) :: spl_local
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: npts_local

        integer, parameter :: niter_old = 3
        integer :: it, ipt
        real(dp) :: t0, t1, dt, sum_y
        real(dp), allocatable :: y_point(:)

        if (npts_local < 1) error stop "bench_old_public: npts_local < 1"

        allocate (y_point(spl_local%num_quantities))

        sum_y = 0.0d0
        best = huge(1.0d0)
        do it = 1, niter_old
            t0 = wall_time()
            do ipt = 1, npts_local
                call evaluate_batch_splines_1d(spl_local, x(ipt), y_point)
                sum_y = sum_y + y_point(1)
            end do
            t1 = wall_time()
            dt = t1 - t0
            best = min(best, dt)
        end do

        print *, "old_public best_s ", best, " pts_per_s ", real(npts_local, dp)/best, &
            " sum_y ", sum_y

        deallocate (y_point)
    end subroutine bench_old_public

    subroutine bench_public(spl_local, x, y, y_expected)
        type(BatchSplineData1D), intent(in) :: spl_local
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout), target :: y(:)
        real(dp), intent(in) :: y_expected(:)

        integer, parameter :: niter_public = 8
        integer :: it
        real(dp) :: t0, t1, dt, warmup
        real(dp) :: diff_max_public

        real(dp), pointer :: y2d(:, :)

        y_ptr = c_loc(y(1))
        call c_f_pointer(y_ptr, y2d, [spl_local%num_quantities, size(x)])

        t0 = wall_time()
        call evaluate_batch_splines_1d_many(spl_local, x, y2d)
        t1 = wall_time()
        warmup = t1 - t0

        best = huge(1.0d0)
        do it = 1, niter_public
            t0 = wall_time()
            call evaluate_batch_splines_1d_many(spl_local, x, y2d)
            t1 = wall_time()
            dt = t1 - t0
            best = min(best, dt)
        end do

        diff_max_public = maxval(abs(y - y_expected))
        print *, "public warmup_s ", warmup
        print *, "public best_s ", best, " pts_per_s ", real(size(x), dp)/best, &
            " max_abs_diff ", diff_max_public
    end subroutine bench_public

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
