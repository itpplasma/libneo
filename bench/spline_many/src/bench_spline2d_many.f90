program bench_spline2d_many
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer
    use batch_interpolate_types, only: BatchSplineData2D
    use batch_interpolate_2d, only: construct_batch_splines_2d, &
                                    construct_batch_splines_2d_lines, &
                                    construct_batch_splines_2d_resident, &
                                    construct_batch_splines_2d_resident_device, &
                                    evaluate_batch_splines_2d, &
                                    evaluate_batch_splines_2d_many, &
                                    destroy_batch_splines_2d
    use spline2d_many_offload, only: spline2d_many_setup, &
                                     spline2d_many_teardown, &
                                     spline2d_many_eval_host, &
                                     spline2d_many_eval_resident
    use util_lcg_rng, only: lcg_state_t, lcg_init, lcg_uniform_0_1
    use util_timer, only: wall_time
    use util_bench_args, only: get_env_int
    implicit none

    integer, parameter :: order(2) = [5, 5]
    integer :: num_points(2)
    integer :: npts
    integer :: niter
    integer :: nbuild
    integer :: nbuild_resident
    integer :: nbuild_repeat
    integer :: neval_repeat
    logical, parameter :: periodic(2) = [.true., .true.]

    real(dp), parameter :: x_min(2) = [1.23d0, -0.7d0]
    real(dp), parameter :: x_max(2) = [x_min(1) + 40.0d0, x_min(2) + 30.0d0]

    type(BatchSplineData2D) :: spl
    real(dp), allocatable, target :: y_ref(:)
    real(dp), allocatable, target :: y_out(:)
    real(dp), allocatable, target :: y_public(:)
    type(c_ptr) :: y_ptr
    real(dp), allocatable :: x_eval(:, :)
    real(dp), allocatable :: y_grid(:, :, :)

    type(lcg_state_t) :: rng
    integer :: num_quantities
    integer :: npts_old
    integer :: i1, i2, iq, ipt
    real(dp) :: best, diff_max

    num_points(1) = get_env_int("LIBNEO_BENCH_N1", 256)
    num_points(2) = get_env_int("LIBNEO_BENCH_N2", 256)
    if (any(num_points < 8)) error stop "bench_spline2d_many: LIBNEO_BENCH_N1/N2 < 8"

    npts = get_env_int("LIBNEO_BENCH_NPTS", 500000)
    if (npts < 1) error stop "bench_spline2d_many: LIBNEO_BENCH_NPTS < 1"

    niter = get_env_int("LIBNEO_BENCH_NITER", 10)
    if (niter < 1) error stop "bench_spline2d_many: LIBNEO_BENCH_NITER < 1"

    nbuild = get_env_int("LIBNEO_BENCH_NBUILD", 2)
    if (nbuild < 1) error stop "bench_spline2d_many: LIBNEO_BENCH_NBUILD < 1"

    nbuild_resident = get_env_int("LIBNEO_BENCH_NBUILD_RESIDENT", nbuild)
    if (nbuild_resident < 1) error stop "bench_spline2d_many: LIBNEO_BENCH_NBUILD_RESIDENT < 1"

    nbuild_repeat = get_env_int("LIBNEO_BENCH_NBUILD_REPEAT", 3)
    if (nbuild_repeat < 1) error stop "bench_spline2d_many: LIBNEO_BENCH_NBUILD_REPEAT < 1"

    neval_repeat = get_env_int("LIBNEO_BENCH_NEVAL_REPEAT", 5)
    if (neval_repeat < 1) error stop "bench_spline2d_many: LIBNEO_BENCH_NEVAL_REPEAT < 1"

    call lcg_init(rng, 24681357_int64)

    num_quantities = get_env_int("LIBNEO_BENCH_NQ", 8)
    if (num_quantities < 1) then
        error stop "bench_spline2d_many: LIBNEO_BENCH_NQ must be >= 1"
    end if

    allocate (y_grid(num_points(1), num_points(2), num_quantities))
    do iq = 1, num_quantities
        do i2 = 1, num_points(2)
            do i1 = 1, num_points(1)
                y_grid(i1, i2, iq) = &
                    cos( &
                    x_min(1) + real(i1 - 1, dp)*(x_max(1) - x_min(1))/ &
                    real(num_points(1) - 1, dp) + real(iq - 1, dp)*0.1d0 &
                    )* &
                    cos( &
                    x_min(2) + real(i2 - 1, dp)*(x_max(2) - x_min(2))/ &
                    real(num_points(2) - 1, dp) + real(iq - 1, dp)*0.2d0 &
                    )
            end do
        end do
    end do

    call bench_build_2d(x_min, x_max, y_grid, spl)
    call bench_build_lines_2d(x_min, x_max, y_grid)
    call bench_build_resident_2d(x_min, x_max, y_grid, spl)
    call bench_build_device_2d()

    allocate (x_eval(2, npts))
    do ipt = 1, npts
        x_eval(1, ipt) = x_min(1) + (x_max(1) - x_min(1))*lcg_uniform_0_1(rng)
        x_eval(2, ipt) = x_min(2) + (x_max(2) - x_min(2))*lcg_uniform_0_1(rng)
    end do

    allocate (y_ref(num_quantities*npts))
    allocate (y_out(num_quantities*npts))
    allocate (y_public(num_quantities*npts))

    call spline2d_many_eval_host(spl%order, spl%num_points, spl%num_quantities, &
                                 spl%periodic, spl%x_min, spl%h_step, spl%coeff, &
                                 x_eval, y_ref)

    print *, "Benchmark 2D parameters"
    print *, "order         ", order
    print *, "num_points    ", num_points
    print *, "num_quantities", num_quantities
    print *, "npts          ", npts
    print *, "niter         ", niter
    print *, "periodic      ", periodic

    npts_old = min(npts, get_env_int("LIBNEO_BENCH_OLD_NPTS", 100000))
    print *, "npts_old      ", npts_old

    call bench_cpu()
    call bench_old_public(npts_old)
#if defined(LIBNEO_ENABLE_OPENACC)
    call bench_openacc()
#endif
    call bench_public()

    call destroy_batch_splines_2d(spl)

contains

    subroutine bench_build_2d(x_min_local, x_max_local, y_grid_local, spl_out)
        real(dp), intent(in) :: x_min_local(2)
        real(dp), intent(in) :: x_max_local(2)
        real(dp), intent(in) :: y_grid_local(:, :, :)
        type(BatchSplineData2D), intent(out) :: spl_out

        integer :: it, rep
        real(dp) :: t0, t1, dt, best_build
        real(dp) :: grid_pts_per_s

        best_build = huge(1.0d0)
        do it = 1, nbuild
            t0 = wall_time()
            do rep = 1, nbuild_repeat
                call construct_batch_splines_2d(x_min_local, x_max_local, &
                                                y_grid_local, order, periodic, spl_out)
                call destroy_batch_splines_2d(spl_out)
            end do
            t1 = wall_time()
            dt = (t1 - t0) / real(nbuild_repeat, dp)
            best_build = min(best_build, dt)
        end do

        call construct_batch_splines_2d(x_min_local, x_max_local, y_grid_local, &
                                        order, periodic, spl_out)

        grid_pts_per_s = real(num_points(1)*num_points(2)*num_quantities, dp) / &
                         best_build
        print *, "build best_s ", best_build, " grid_pts_per_s ", grid_pts_per_s
    end subroutine bench_build_2d

    subroutine bench_build_lines_2d(x_min_local, x_max_local, y_grid_local)
        real(dp), intent(in) :: x_min_local(2)
        real(dp), intent(in) :: x_max_local(2)
        real(dp), intent(in) :: y_grid_local(:, :, :)

        type(BatchSplineData2D) :: spl_tmp
        integer :: it, rep
        real(dp) :: t0, t1, dt, best_build
        real(dp) :: grid_pts_per_s

        best_build = huge(1.0d0)
        do it = 1, nbuild
            t0 = wall_time()
            do rep = 1, nbuild_repeat
                call construct_batch_splines_2d_lines(x_min_local, x_max_local, &
                                                      y_grid_local, order, periodic, &
                                                      spl_tmp)
                call destroy_batch_splines_2d(spl_tmp)
            end do
            t1 = wall_time()
            dt = (t1 - t0) / real(nbuild_repeat, dp)
            best_build = min(best_build, dt)
        end do

        grid_pts_per_s = real(num_points(1)*num_points(2)*num_quantities, dp) / &
                         best_build
        print *, "build_lines best_s ", best_build, " grid_pts_per_s ", grid_pts_per_s
    end subroutine bench_build_lines_2d

    subroutine bench_build_resident_2d(x_min_local, x_max_local, y_grid_local, spl_out)
        real(dp), intent(in) :: x_min_local(2)
        real(dp), intent(in) :: x_max_local(2)
        real(dp), intent(in) :: y_grid_local(:, :, :)
        type(BatchSplineData2D), intent(inout) :: spl_out

        integer :: it, rep
        real(dp) :: t0, t1, dt, best_build
        real(dp) :: grid_pts_per_s

#if defined(LIBNEO_ENABLE_OPENACC)
        call destroy_batch_splines_2d(spl_out)

        best_build = huge(1.0d0)
        do it = 1, nbuild_resident
            t0 = wall_time()
            do rep = 1, nbuild_repeat
                call construct_batch_splines_2d_resident(x_min_local, x_max_local, &
                                                         y_grid_local, &
                                                         order, periodic, spl_out)
                !$acc wait
                call destroy_batch_splines_2d(spl_out)
            end do
            t1 = wall_time()
            dt = (t1 - t0) / real(nbuild_repeat, dp)
            best_build = min(best_build, dt)
        end do

        call construct_batch_splines_2d_resident(x_min_local, x_max_local, &
                                                 y_grid_local, order, periodic, spl_out)

        grid_pts_per_s = real(num_points(1)*num_points(2)*num_quantities, dp) / &
                         best_build
        print *, "build_resident best_s ", best_build, " grid_pts_per_s ", &
            grid_pts_per_s
#else
        call destroy_batch_splines_2d(spl_out)
        call construct_batch_splines_2d(x_min_local, x_max_local, y_grid_local, &
                                        order, periodic, spl_out)
#endif
    end subroutine bench_build_resident_2d

    subroutine bench_build_device_2d()
        type(BatchSplineData2D), save :: spl_dev
        integer :: it, rep
        real(dp) :: t0, t1, dt, best_build
        real(dp) :: grid_pts_per_s

#if defined(LIBNEO_ENABLE_OPENACC)
        best_build = huge(1.0d0)
        do it = 1, nbuild_resident
            t0 = wall_time()
            do rep = 1, nbuild_repeat
                call construct_batch_splines_2d_resident_device(&
                    x_min, x_max, y_grid, order, periodic, spl_dev, &
                    update_host=.true., assume_y_present=.false.)
                !$acc wait
                call destroy_batch_splines_2d(spl_dev)
            end do
            t1 = wall_time()
            dt = (t1 - t0) / real(nbuild_repeat, dp)
            best_build = min(best_build, dt)
        end do
        grid_pts_per_s = real(num_points(1)*num_points(2)*num_quantities, dp) / &
                         best_build
        print *, "build_device_host best_s ", best_build, " grid_pts_per_s ", &
            grid_pts_per_s

        !$acc enter data copyin(y_grid)
        best_build = huge(1.0d0)
        do it = 1, nbuild_resident
            t0 = wall_time()
            do rep = 1, nbuild_repeat
                call construct_batch_splines_2d_resident_device(&
                    x_min, x_max, y_grid, order, periodic, spl_dev, &
                    update_host=.false., assume_y_present=.true.)
                !$acc wait
                call destroy_batch_splines_2d(spl_dev)
            end do
            t1 = wall_time()
            dt = (t1 - t0) / real(nbuild_repeat, dp)
            best_build = min(best_build, dt)
        end do
        !$acc exit data delete(y_grid)
        grid_pts_per_s = real(num_points(1)*num_points(2)*num_quantities, dp) / &
                         best_build
        print *, "build_device_gpu best_s ", best_build, " grid_pts_per_s ", &
            grid_pts_per_s
#endif
    end subroutine bench_build_device_2d

    subroutine bench_cpu()
        integer :: it, rep
        real(dp) :: t0, t1, dt

        best = huge(1.0d0)
        do it = 1, niter
            t0 = wall_time()
            do rep = 1, neval_repeat
                call spline2d_many_eval_host(spl%order, spl%num_points, &
                                             spl%num_quantities, spl%periodic, &
                                             spl%x_min, spl%h_step, spl%coeff, &
                                             x_eval, y_out)
            end do
            t1 = wall_time()
            dt = (t1 - t0) / real(neval_repeat, dp)
            best = min(best, dt)
        end do

        diff_max = maxval(abs(y_out - y_ref))
        print *, "cpu best_s ", best, " pts_per_s ", real(npts, dp)/best, &
            " max_abs_diff ", &
            diff_max
    end subroutine bench_cpu

    subroutine bench_old_public(npts_local)
        integer, intent(in) :: npts_local

        integer, parameter :: niter_old = 2
        integer :: it, ipt_local
        real(dp) :: t0, t1, dt, sum_y
        real(dp), allocatable :: y_point(:)

        if (npts_local < 1) error stop "bench_old_public: npts_local < 1"
        allocate (y_point(spl%num_quantities))

        sum_y = 0.0d0
        best = huge(1.0d0)
        do it = 1, niter_old
            t0 = wall_time()
            do ipt_local = 1, npts_local
                call evaluate_batch_splines_2d(spl, x_eval(:, ipt_local), y_point)
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

    subroutine bench_public()
        integer, parameter :: niter_public = 6
        integer :: it
        real(dp) :: t0, t1, dt, warmup
        real(dp) :: diff_max_public

        real(dp), pointer :: y2d(:, :)

        y_ptr = c_loc(y_public(1))
        call c_f_pointer(y_ptr, y2d, [spl%num_quantities, npts])

        t0 = wall_time()
        call evaluate_batch_splines_2d_many(spl, x_eval, y2d)
        t1 = wall_time()
        warmup = t1 - t0

        best = huge(1.0d0)
        do it = 1, niter_public
            t0 = wall_time()
            call evaluate_batch_splines_2d_many(spl, x_eval, y2d)
            t1 = wall_time()
            dt = t1 - t0
            best = min(best, dt)
        end do

        diff_max_public = maxval(abs(y_public - y_ref))
        print *, "public warmup_s ", warmup
        print *, "public best_s ", best, " pts_per_s ", real(npts, dp)/best, &
            " max_abs_diff ", diff_max_public
    end subroutine bench_public

    subroutine bench_openacc()
        integer :: it, rep
        real(dp) :: t0, t1, dt
        real(dp) :: tsetup0, tsetup1

        tsetup0 = wall_time()
        call spline2d_many_setup(spl%coeff, x_eval, y_out)
        !$acc wait
        tsetup1 = wall_time()
        print *, "openacc setup_s ", tsetup1 - tsetup0
        call spline2d_many_eval_resident(spl%order, spl%num_points, &
                                         spl%num_quantities, spl%periodic, &
                                         spl%x_min, spl%h_step, spl%coeff, &
                                         x_eval, y_out)
        !$acc wait

        best = huge(1.0d0)
        do it = 1, niter
            t0 = wall_time()
            do rep = 1, neval_repeat
                call spline2d_many_eval_resident(spl%order, spl%num_points, &
                                                 spl%num_quantities, spl%periodic, &
                                                 spl%x_min, spl%h_step, spl%coeff, &
                                                 x_eval, y_out)
            end do
            !$acc wait
            t1 = wall_time()
            dt = (t1 - t0) / real(neval_repeat, dp)
            best = min(best, dt)
        end do

        !$acc update self(y_out)
        diff_max = maxval(abs(y_out - y_ref))
        print *, "openacc best_s ", best, " pts_per_s ", real(npts, dp)/best, &
            " max_abs_diff ", diff_max
        call spline2d_many_teardown(spl%coeff, x_eval, y_out)
    end subroutine bench_openacc

end program bench_spline2d_many
